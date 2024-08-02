# source: DART-eval


rule dsqtl_process:
    input:
        "results/dsqtl/yoruban.dsqtls.benchmarking.tsv.gz",
        "results/genome.fa.gz",
    output:
        "results/dsqtl/processed/test.parquet",
    run:
        V = (
            pd.read_csv(
                input[0], sep="\t",
                usecols=[
                    "var.chrom", "var.pos", "var.allele1", "var.allele2", "var.label"
                ]
            )
            .rename(columns={
                "var.chrom": "chrom", "var.pos": "pos",
                "var.allele1": "ref", "var.allele2": "alt", "var.label": "label"
            })
        )
        V.chrom = V.chrom.str.replace("chr", "")
        V["label"] = V.label == 1
        print(V.shape)
        V = lift_hg19_to_hg38(V)
        V = V[V.pos != -1]
        print(V.shape)
        V = V[(V.ref.str.len()==1) & (V.alt.str.len()==1)]
        print(V.shape)
        genome = Genome(input[1])
        V = check_ref_alt(V, genome)
        print(V.shape)
        V = sort_chrom_pos(V)
        V.to_parquet(output[0], index=False) 
