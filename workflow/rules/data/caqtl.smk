# source: DART-eval


rule caqtl_process:
    input:
        "results/caqtl/Afr.CaQTLS.tsv.gz",
        "results/genome.fa.gz",
    output:
        "results/caqtl/processed/test.parquet",
    run:
        V = (
            pd.read_csv(
                input[0], sep="\t",
                usecols=[
                    "chr_hg38", "pos_hg38", "allele1", "allele2", "label"
                ]
            )
            .rename(columns={
                "chr_hg38": "chrom", "pos_hg38": "pos",
                "allele1": "ref", "allele2": "alt",
            })
        )
        V.chrom = V.chrom.str.replace("chr", "")
        print(V.shape)
        V = V[(V.ref.str.len()==1) & (V.alt.str.len()==1)]
        print(V.shape)
        genome = Genome(input[1])
        V = check_ref_alt(V, genome)
        print(V.shape)
        V.to_parquet(output[0], index=False) 


rule caqtl_subsample:
    input:
        "results/caqtl/processed/test.parquet",
    output:
        "results/caqtl/subsampled/test.parquet",
    run:
        df = pd.read_parquet(input[0])
        all_pos = []
        all_neg_matched = []
        for chrom in tqdm(df.chrom.unique()):
            df_c = df[df.chrom == chrom]
            pos = df_c[df_c.label]
            neg = df_c[~df_c.label]
            neg = neg.sample(len(pos), random_state=42)
            all_pos.append(pos)
            all_neg_matched.append(neg)
    
        pos = pd.concat(all_pos, ignore_index=True)
        pos["match_group"] = np.arange(len(pos))
        neg_matched = pd.concat(all_neg_matched, ignore_index=True)
        neg_matched["match_group"] = np.arange(len(neg_matched))
        res = pd.concat([pos, neg_matched], ignore_index=True)
        res = sort_chrom_pos(res)
        res.to_parquet(output[0], index=False)
