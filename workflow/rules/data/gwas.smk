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


#rule gwas_filt:
#    input:
#        "results/gwas/processed.parquet",
#    output:
#        "results/gwas/filt.parquet",
#    run:
#        V = pd.read_parquet(input[0])
#        V = V[V.method=="SUSIE"]
#        V = V[(~V.LD_HWE) & (~V.LD_SV)]
#        # reinterpreting trait as "causal in trait" rather than "tested" in trait
#        V.loc[V.pip <= 0.9, "trait"] = ""
#        V = V.groupby(COORDINATES).agg({
#            "pip": "max", "maf": "mean", "trait": "unique",
#        }).reset_index()
#        V.loc[V.pip > 0.9, "label"] = True
#        V.loc[V.pip < 0.01, "label"] = False
#        V = V.dropna(subset="label")
#        V.trait = V.trait.progress_apply(
#            lambda traits: ",".join(sorted([trait for trait in traits if trait != ""]))
#        )
#        print(V)
#        V.to_parquet(output[0], index=False)
#
#
#rule gwas_intermediate:
#    input:
#        "results/gwas/processed.parquet",
#    output:
#        "results/gwas/intermediate/test.parquet",
#    run:
#        V = pd.read_parquet(input[0])
#        V = V[V.method=="SUSIE"]
#        V = V[(~V.LD_HWE) & (~V.LD_SV)]
#        V = V.groupby(COORDINATES).agg({"pip": "max"}).reset_index()
#        V = V[V.pip.between(0.01, 0.9)]
#        print(V)
#        V.to_parquet(output[0], index=False)
#