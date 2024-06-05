rule splicing_assay_process:
    input:
        "results/data/splicing_assay/Pangolin_pred.csv",
        "results/genome.fa.gz",
    output:
        "results/data/splicing_assay/test.parquet",
        "results/features/gonzalobenegas/splicing_assay/Pangolin.parquet",
    run:
        V = pd.read_csv(input[0]).rename(columns={
            "CHROM": "chrom",
            "POS": "pos",
            "REF": "ref",
            "ALT": "alt",
            "strong_lof": "label",
            "pan_diff": "Pangolin",
        })
        print(V)
        V.chrom = V.chrom.str.replace("chr", "")
        V = filter_snp(V)
        V = lift_hg19_to_hg38(V)
        V = V[V.pos != -1]
        genome = Genome(input[1])
        V = check_ref_alt(V, genome)
        V = sort_chrom_pos_ref_alt(V)
        print(V)
        V[COORDINATES + ["label", "gene", "exon/intron", "min_distance"]].to_parquet(output[0], index=False)
        V[["Pangolin"]].to_parquet(output[1], index=False)
