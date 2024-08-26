# all hacky for now


rule ldscore_process:
    input:
        "workflow/notebooks/UKBB.ALL.ldscore/UKBB.EUR.l2.ldscore.gz",
        "results/genome.fa.gz",
    output:
        "results/ldscore/processed.parquet",
    run:
        V = (
            pl.read_csv(
                input[0], separator="\t", columns=["SNP", "L2"],
            )
            .with_columns(
                pl.col("SNP").str.split_exact(":", 3)
                .struct.rename_fields(COORDINATES)
            )
            .with_columns(
                pl.col("SNP").struct.field("chrom"),
                pl.col("SNP").struct.field("pos").cast(int),
                pl.col("SNP").struct.field("ref"),
                pl.col("SNP").struct.field("alt"),
            )
            .drop("SNP")
            .select(COORDINATES + ["L2"])
            .to_pandas()
        )
        print(V)
        V = filter_snp(V)
        print(V.shape)
        V = lift_hg19_to_hg38(V)
        V = V[V.pos != -1]
        print(V.shape)
        genome = Genome(input[1])
        V = check_ref_alt(V, genome)
        print(V.shape)
        V = sort_variants(V)
        V.to_parquet(output[0], index=False)


rule ldscore_feature:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/ldscore/processed.parquet",
    output:
        "results/dataset/{dataset}/features/LDScore.parquet",
    run:
        V = pd.read_parquet(input[0])
        ldscore = pd.read_parquet(input[1]).rename(columns={"L2": "score"})
        V = V.merge(ldscore, on=COORDINATES, how="left")
        V[["score"]].to_parquet(output[0], index=False)
