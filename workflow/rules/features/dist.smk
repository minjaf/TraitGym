rule tss_dist_feature:
    input:
        "results/tss.parquet",
    output:
        "results/features/{dataset}/minus_TSS_dist.parquet",
    run:
        V = load_dataset(wildcards["dataset"], split="test").to_pandas()
        V["start"] = V["pos"] - 1
        V["end"] = V["pos"]
        tss = pd.read_parquet(input[0], columns=["chrom", "start", "end"])
        V = bf.closest(V, tss)
        V["minus_tss_dist"] = -V.distance
        V = sort_chrom_pos(V)
        V[["minus_tss_dist"]].to_parquet(output[0], index=False)
