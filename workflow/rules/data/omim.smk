# Curated regulatory OMIM variants, Table S6 from:
# Smedley, Damian, et al. "A whole-genome analysis framework for effective
# identification of pathogenic regulatory variants in Mendelian disease." The
# American Journal of Human Genetics 99.3 (2016): 595-606.


rule omim_download:
    output:
        temp("results/omim/variants.xslx"),
    shell:
        "wget -O {output} https://ars.els-cdn.com/content/image/1-s2.0-S0002929716302786-mmc2.xlsx"


rule omim_process:
    input:
        "results/omim/variants.xslx",
    output:
        "results/omim/variants.parquet",
    run:
        xls = pd.ExcelFile(input[0])
        sheet_names = xls.sheet_names

        dfs = []
        for variant_type in sheet_names[1:]:
            df = pd.read_excel(input[0], sheet_name=variant_type)
            dfs.append(df)
        V = pd.concat(dfs)
        V = V[["Chr", "Position", "Ref", "Alt", "OMIM"]].rename(columns={
            "Chr": "chrom", "Position": "pos", "Ref": "ref", "Alt": "alt"
        })
        V.chrom = V.chrom.str.replace("chr", "")
        V = filter_snp(V)
        V = lift_hg19_to_hg38(V)
        V = V[V.pos != -1]
        print(V)
        V.to_parquet(output[0], index=False)


rule omim_dataset:
    input:
        "results/omim/variants.annot_with_cre.parquet",
        "results/gnomad/common.parquet",
    output:
        "results/dataset/omim_subsampled_{k,\d+}/test.parquet",
    run:
        k = int(wildcards.k)
        pos = pd.read_parquet(input[0])
        pos["label"] = True
        neg = pd.read_parquet(input[1])
        neg["label"] = False
        all_pos = []
        all_neg = []
        for c in tqdm(pos.consequence.unique()):
            for chroms in ODD_EVEN_CHROMS:
                pos_c = pos[(pos.consequence == c) & (pos.chrom.isin(chroms))]
                neg_c = neg[(neg.consequence == c) & (neg.chrom.isin(chroms))]
                if len(neg_c) < len(pos_c) * k:
                    print(f"Subsampling {c}")
                    pos_c = pos_c.sample(n=len(neg_c) // k, random_state=42)
                all_pos.append(pos_c)
                all_neg.append(neg_c.sample(n=len(pos_c) * k, random_state=42))
        pos = pd.concat(all_pos, ignore_index=True)
        pos["match_group"] = np.arange(len(pos))
        neg = pd.concat(all_neg, ignore_index=True)
        neg["match_group"] = np.repeat(pos.match_group.values, k)
        V = pd.concat([pos, neg], ignore_index=True)
        V = sort_variants(V)
        print(V)
        V.to_parquet(output[0], index=False)
