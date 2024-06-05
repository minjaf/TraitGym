rule dataset_to_vcf:
    output:
        "results/data/{dataset}.vcf.gz",
    run:
        V = load_dataset(wildcards.dataset, split="test").to_pandas()
        V["id"] = "."
        V[["chrom", "pos", "id", "ref", "alt"]].to_csv(
            output[0], sep="\t", index=False, header=False,
        )


# obtained from the website
rule cadd_process:
    input:
        "results/data/cadd/{dataset}.tsv.gz",
    output:
        "results/features/{dataset}/CADD.parquet",
    run:
        V = load_dataset(wildcards.dataset, split="test").to_pandas()
        scores = pd.read_csv(
            input[0], sep="\t", comment="#", header=None,
            names=["chrom", "pos", "ref", "alt", "RawScore", "PHRED"],
            dtype={"chrom": "str"},
        )
        # make sure all on the same half-plane
        scores.RawScore = scores.RawScore - scores.RawScore.min()
        V = V.merge(scores, how="left", on=COORDINATES)
        print(V)
        V[["RawScore"]].to_parquet(output[0], index=False)
