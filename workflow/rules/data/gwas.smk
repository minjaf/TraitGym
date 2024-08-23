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


rule gwas_process:
    input:
        "results/gwas/raw/release1.1/UKBB_94traits_release1.bed.gz",
        "results/gwas/raw/release1.1/UKBB_94traits_release1_regions.bed.gz",
        "results/genome.fa.gz",
    output:
        "results/gwas/processed.parquet",
    run:
        V1 = (
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
            .group_by(["chrom", "pos", "ref", "alt", "trait"])
            .agg(
                pl.mean("pip"),
                (pl.max("pip") - pl.min("pip")).alias("pip_diff"),
                pl.count().alias("pip_n"),
            )
            .filter(pl.col("pip_n") == 2, pl.col("pip_diff") < 0.05)
            .drop("pip_n", "pip_diff")
            # TODO: should not drop until we have the variants with pip around 0
            # in the second file. as not to drop such control variants
            # maybe need to do two copies of the second file, one for susie and one for
            # finemap
        )
        print(V1)
        print(V1.filter(pl.col("pip") > 0.9).unique(COORDINATES))
        raise Exception("debug")
        V.chrom = V.chrom.str.replace("chr", "")
        print(V.shape)
        V = lift_hg19_to_hg38(V)
        V = V[V.pos != -1]
        print(V.shape)
        V = V[(V.ref.str.len()==1) & (V.alt.str.len()==1)]
        print(V.shape)
        genome = Genome(input[4])
        V = check_ref_alt(V, genome)
        print(V.shape)
        V.to_parquet(output[0], index=False) 
#
#
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