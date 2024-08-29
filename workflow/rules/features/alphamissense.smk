rule alphamissense_download:
    output:
        temp("results/alphamissense/scores.tsv.gz"),
    shell:
        "wget -O {output} https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz"


rule alphamissense_process:
    input:
        "results/alphamissense/scores.tsv.gz",
    output:
        "results/alphamissense/scores.parquet",
    run:
        V = (
            pl.read_csv(
                input[0], separator="\t", skip_rows=3,
                columns=["#CHROM", "POS", "REF", "ALT", "am_pathogenicity"],
            )
            .rename({
                "#CHROM": "chrom",
                "POS": "pos",
                "REF": "ref",
                "ALT": "alt",
                "am_pathogenicity": "score",
            })
            .with_columns(pl.col("chrom").str.replace("chr", ""))
            # consequence in multiple isoforms: take the most pathogenic
            .group_by(COORDINATES).agg(pl.max("score"))
            .sort(COORDINATES)
        )
        print(V)
        V.write_parquet(output[0])


rule alphamissense_features:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/alphamissense/scores.parquet",
    output:
        "results/dataset/{dataset}/features/AlphaMissense.parquet",
    run:
        V = pd.read_parquet(input[0], columns=COORDINATES)
        score = pd.read_parquet(input[1])
        V = V.merge(score, how="left", on=COORDINATES)
        print(V)
        V[["score"]].to_parquet(output[0], index=False)
