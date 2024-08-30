rule gwas_download:
    output:
        temp("results/gwas/UKBB_94traits_release1.1.tar.gz"),
        "results/gwas/raw/release1.1/UKBB_94traits_release1.bed.gz",
        "results/gwas/raw/release1.1/UKBB_94traits_release1.cols",
        "results/gwas/raw/release1.1/UKBB_94traits_release1_regions.bed.gz",
        "results/gwas/raw/release1.1/UKBB_94traits_release1_regions.cols",
    params:
        directory("results/gwas/raw"),
    shell:
        """
        wget -O {output[0]} https://www.dropbox.com/s/cdsdgwxkxkcq8cn/UKBB_94traits_release1.1.tar.gz?dl=1 &&
        mkdir -p {params} &&
        tar -xzvf {output[0]} -C {params}
        """


rule gwas_process_main_file:
    input:
        "results/gwas/raw/release1.1/UKBB_94traits_release1.bed.gz",
    output:
        "results/gwas/main_file.parquet",
    run:
        V = (
            pl.read_csv(
                input[0], separator="\t", has_header=False,
                columns=[0, 2, 5, 6, 10, 11, 17, 21, 22],
                new_columns=[
                    "chrom", "pos", "ref", "alt", "method", "trait", "pip", "LD_HWE",
                    "LD_SV",
                ],
                schema_overrides={"column_3": float}
            )
            .with_columns(pl.col("pos").cast(int))
            .filter(~pl.col("LD_HWE"), ~pl.col("LD_SV"))
            .select(["chrom", "pos", "ref", "alt", "trait", "method", "pip"])
        )
        print(V)
        V.write_parquet(output[0])


rule gwas_process_secondary_file:
    input:
        "results/gwas/raw/release1.1/UKBB_94traits_release1_regions.bed.gz",
    output:
        "results/gwas/secondary_file.parquet",
    run:
        V = (
            pl.read_csv(
                input[0], separator="\t", has_header=False, columns=[4, 6, 7, 8],
                new_columns=["trait", "variant", "success_finemap", "success_susie"],
            )
            .filter(pl.col("success_finemap"), pl.col("success_susie"))
            .drop("success_finemap", "success_susie")
            .with_columns(
                pl.col("variant").str.split_exact(":", 3)
                .struct.rename_fields(COORDINATES)
            )
            .with_columns(
                pl.col("variant").struct.field("chrom"),
                pl.col("variant").struct.field("pos").cast(int),
                pl.col("variant").struct.field("ref"),
                pl.col("variant").struct.field("alt"),
            )
            .drop("variant")
            .select(["chrom", "pos", "ref", "alt", "trait"])
            .with_columns(pl.lit(0.0).alias("pip"))
        )
        print(V)
        V = (
            pl.concat([
                V.with_columns(pl.lit("SUSIE").alias("method")),
                V.with_columns(pl.lit("FINEMAP").alias("method")),
            ])
            .select(["chrom", "pos", "ref", "alt", "trait", "method", "pip"])
        )
        print(V)
        V.write_parquet(output[0])


rule gwas_process:
    input:
        "results/gwas/main_file.parquet",
        "results/gwas/secondary_file.parquet",
        "results/genome.fa.gz",
    output:
        "results/gwas/processed.parquet",
    run:
        V = (
            pl.concat([pl.read_parquet(input[0]), pl.read_parquet(input[1])])
            .unique(
                ["chrom", "pos", "ref", "alt", "trait", "method"],
                keep="first", maintain_order=True
            )
            .group_by(["chrom", "pos", "ref", "alt", "trait"])
            .agg(
                pl.mean("pip"),
                (pl.max("pip") - pl.min("pip")).alias("pip_diff"),
                pl.count().alias("pip_n"),
            )
            .filter(pl.col("pip_n") == 2, pl.col("pip_diff") < 0.05)
            .drop("pip_n", "pip_diff")
        )
        print(V)
        V = (
            V.with_columns(
                pl.when(pl.col("pip") > 0.9)
                .then(pl.col("trait"))
                .otherwise(pl.lit(None))
                .alias("trait")
            )
            .group_by(COORDINATES)
            .agg(pl.max("pip"), pl.col("trait").drop_nulls().unique())
            .with_columns(pl.col("trait").list.sort().list.join(","))
            .to_pandas()
        )
        print(V)
        V.chrom = V.chrom.str.replace("chr", "")
        V = filter_snp(V)
        print(V.shape)
        V = lift_hg19_to_hg38(V)
        V = V[V.pos != -1]
        print(V.shape)
        genome = Genome(input[2])
        V = check_ref_alt(V, genome)
        print(V.shape)
        V = sort_variants(V)
        print(V)
        V.to_parquet(output[0], index=False) 


rule gwas_match:
    input:
        "results/gwas/processed.parquet",
        "results/ldscore/UKBB.EUR.ldscore.annot_with_cre.parquet",
        "results/tss.parquet",
    output:
        "results/dataset/gwas_matched_{k,\d+}/test.parquet",
    run:
        k = int(wildcards.k)
        V = (
            pl.read_parquet(input[0])
            .with_columns(
                pl.when(pl.col("pip") > 0.9).then(True)
                .when(pl.col("pip") < 0.01).then(False)
                .otherwise(None)
                .alias("label")
            )
            .drop_nulls()
            .to_pandas()
        )

        annot = pd.read_parquet(input[1])
        V = V.merge(annot, on=COORDINATES, how="inner")

        V["start"] = V.pos - 1
        V["end"] = V.pos

        tss = pd.read_parquet(input[2], columns=["chrom", "start", "end"])

        V = bf.closest(V, tss).rename(columns={
            "distance": "tss_dist"
        }).drop(columns=["start", "end", "chrom_", "start_", "end_"])

        match_features = ["maf", "ld_score", "tss_dist"]

        consequences = V[V.label].consequence.unique()
        V_cs = []
        for c in consequences:
            print(c)
            V_c = V[V.consequence == c].copy()
            for f in match_features:
                V_c[f"{f}_scaled"] = RobustScaler().fit_transform(V_c[f].values.reshape(-1, 1))
            print(V_c.label.value_counts())
            V_c = match_columns_k(V_c, "label", [f"{f}_scaled" for f in match_features], k)
            V_c["match_group"] = c + "_" + V_c.match_group.astype(str)
            print(V_c.label.value_counts())
            print(V_c.groupby("label")[match_features].median())
            V_c.drop(columns=[f"{f}_scaled" for f in match_features], inplace=True)
            V_cs.append(V_c)
        V = pd.concat(V_cs, ignore_index=True)
        V = sort_variants(V)
        print(V)
        V.to_parquet(output[0], index=False)
