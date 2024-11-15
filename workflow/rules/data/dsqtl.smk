# source: DART-eval

# TODO: check if this can be downloaded from the original paper (I think
# the SVM paper)

# Supplementary Table 1
# Lee, Dongwon, et al. "A method to predict the impact of regulatory variants
# from DNA sequence." Nature Genetics 47.8 (2015): 955-961.

# https://static-content.springer.com/esm/art%3A10.1038%2Fng.3331/MediaObjects/41588_2015_BFng3331_MOESM26_ESM.xlsx


rule dsqtl_download:
    output:
        "results/dsqtl/variants.xlsx",
    shell:
        "wget -O {output} https://static-content.springer.com/esm/art%3A10.1038%2Fng.3331/MediaObjects/41588_2015_BFng3331_MOESM26_ESM.xlsx"


rule dsqtl_process:
    input:
        "results/dsqtl/variants.xlsx",
        "results/genome.fa.gz",
    output:
        "results/dsqtl/variants.parquet",
    run:
        V = (
            pd.read_excel(
                input[0], sheet_name="SuppTable1",
                usecols=[
                    "chrom_hg19", "pos_hg19", "allele1", "allele2", "label",
                ]
            )
            .rename(columns={
                "chrom_hg19": "chrom", "pos_hg19": "pos",
                "allele1": "ref", "allele2": "alt",
            })
        )
        print(V)
        V.chrom = V.chrom.str.replace("chr", "")
        V["label"] = V.label == 1
        V = filter_snp(V)
        print(V.shape)
        V = lift_hg19_to_hg38(V)
        V = V[V.pos != -1]
        print(V.shape)
        genome = Genome(input[1])
        V = check_ref_alt(V, genome)
        print(V.shape)
        V = sort_variants(V)
        print(V)
        V.to_parquet(output[0], index=False) 


rule dsqtl_dataset_intermediate:
    input:
        "results/dsqtl/variants.annot.parquet",
        "results/tss.parquet",
    output:
        "results/intermediate_dataset/dsqtl_proximal.parquet",
        "results/intermediate_dataset/dsqtl_distal.parquet",
    run:
        V = pd.read_parquet(input[0])
        V = V.drop_duplicates(COORDINATES)
        V = V[V.consequence.isin(NON_EXONIC_FULL)]
        assert len(V) == len(V.drop_duplicates(COORDINATES))
        V["start"] = V.pos - 1
        V["end"] = V.pos
        tss = pd.read_parquet(input[1])
        V = bf.closest(V, tss).rename(columns={
            "distance": "tss_dist", "gene_id_": "gene",
        }).drop(columns=["start", "end", "chrom_", "start_", "end_"])
        print(V.tss_dist.dtype)
        assert V.tss_dist.notna().all()
        V.tss_dist = V.tss_dist.astype(int)
        print(V.tss_dist.dtype)
        V["proximal"] = V.tss_dist < 1000
        V[V.proximal].drop(columns="proximal").to_parquet(output[0], index=False)
        V[~V.proximal].drop(columns="proximal").to_parquet(output[1], index=False)


rule dsqtl_dataset:
    input:
        "results/intermediate_dataset/dsqtl_{proximity}.parquet",
    output:
        "results/dataset/dsqtl_{proximity}_matched_{k,\d+}/test.parquet",
    run:
        k = int(wildcards.k)
        V = pd.read_parquet(input[0])
        V["super_proximal"] = V.tss_dist < 100
        cols = [
            "gene", "super_proximal",
        ]
        print(V)
        V = match_cols(V[V.label], V[~V.label], cols, k=k, minimize_dist_col="tss_dist")
        V = V.drop(columns=["super_proximal"])
        print(V)
        print(V.label.sum())
        print(
            V.label.mean(),
            average_precision_score(V.label, -V.tss_dist),
        )
        V.to_parquet(output[0], index=False)