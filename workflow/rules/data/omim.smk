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


rule mendelian_traits_dataset_intermediate:
    input:
        "results/omim/variants.annot_with_cre.annot_MAF.parquet",
        "results/gnomad/MAF_above_0.1.annot_with_cre.parquet",
        "results/tss.parquet",
    output:
        "results/intermediate_dataset/mendelian_traits_proximal.parquet",
        "results/intermediate_dataset/mendelian_traits_distal.parquet",
    run:
        pos = pd.read_parquet(input[0])
        pos.maf = pos.maf.fillna(0)
        pos = pos[pos.maf < 1 / 100]
        pos["label"] = True
        neg = pd.read_parquet(input[1], columns=COORDINATES + ["consequence", "MAF"])
        neg = neg.rename(columns={"MAF": "maf"})
        neg = neg[neg.maf > 1 / 100]
        neg["label"] = False
        V = pd.concat([pos, neg], ignore_index=True)
        V = V[V.consequence.isin(NON_EXONIC_FULL)]
        assert len(V) == len(V.drop_duplicates(COORDINATES))
        V["start"] = V.pos - 1
        V["end"] = V.pos
        tss = pd.read_parquet(input[2])
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


#rule mendelian_traits_dataset:
#    input:
#        "results/intermediate_dataset/mendelian_traits_{proximity}.parquet",
#    output:
#        "results/dataset/mendelian_traits_{proximity,proximal|distal}_matched_{k,\d+}/test.parquet",
#    run:
#        k = int(wildcards.k)
#        V = pd.read_parquet(input[0])
#        V["super_proximal"] = V.tss_dist < 100
#        cols = ["gene", "super_proximal"]
#        V = match_cols(V[V.label], V[~V.label], cols, k=k, minimize_dist_col="tss_dist")
#        V = V.drop(columns=["super_proximal"])
#        print(V)
#        print(V.label.sum())
#        print(V.label.mean(), average_precision_score(V.label, -V.tss_dist))
#        V.to_parquet(output[0], index=False)


#rule mendelian_traits_subset_trait:
#    input:
#        "results/dataset/mendelian_traits_matched_{k}/test.parquet",
#    output:
#        "results/dataset/mendelian_traits_matched_{k}/subset/{t}.parquet",
#    wildcard_constraints:
#        t="|".join(select_omim_traits),
#    run:
#        V = pd.read_parquet(input[0])
#        target_size = 1 + int(wildcards.k)
#        V = V[(~V.label) | (V.OMIM == f"MIM {wildcards.t}")]
#        match_group_size = V.match_group.value_counts() 
#        match_groups = match_group_size[match_group_size == target_size].index
#        V = V[V.match_group.isin(match_groups)]
#        V[COORDINATES].to_parquet(output[0], index=False)


rule mendelian_traits_all_intermediate_dataset:
    input:
        "results/omim/variants.annot_with_cre.annot_MAF.parquet",
        "results/gnomad/MAF_above_0.1.annot_with_cre.parquet",
        "results/tss.parquet",
    output:
        "results/intermediate_dataset/mendelian_traits_proximal_all.parquet",
        "results/intermediate_dataset/mendelian_traits_distal_all.parquet",
    run:
        pos = pd.read_parquet(input[0])
        pos.maf = pos.maf.fillna(0)
        pos = pos[pos.maf < 1 / 100]
        pos = pos[pos.consequence.isin(NON_EXONIC_FULL)]
        pos["label"] = True
        neg = pd.read_parquet(input[1], columns=COORDINATES + ["consequence", "MAF"])
        neg = neg.rename(columns={"MAF": "maf"})
        neg = neg[neg.maf > 5 / 100]
        neg = neg[neg.chrom.isin(pos.chrom.unique())]
        neg = neg[neg.consequence.isin(NON_EXONIC_FULL)]
        neg["label"] = False
        V = pd.concat([pos, neg], ignore_index=True)
        assert len(V) == len(V.drop_duplicates(COORDINATES))
        V["start"] = V.pos - 1
        V["end"] = V.pos
        tss = pd.read_parquet(input[2])
        V = bf.closest(V, tss).rename(columns={
            "distance": "tss_dist", "gene_id_": "gene",
        }).drop(columns=["start", "end", "chrom_", "start_", "end_"])
        print(V.tss_dist.dtype)
        assert V.tss_dist.notna().all()
        V.tss_dist = V.tss_dist.astype(int)
        print(V.tss_dist.dtype)
        V = sort_variants(V)
        V["proximal"] = V.tss_dist < 1000
        V[V.proximal].drop(columns="proximal").to_parquet(output[0], index=False)
        V[~V.proximal].drop(columns="proximal").to_parquet(output[1], index=False)


#rule mendelian_traits_all_dataset:
#    input:
#        "results/intermediate_dataset/mendelian_traits_{proximity}_all.parquet",
#    output:
#        "results/dataset/mendelian_traits_{proximity,proximal|distal}_all_matched_{k,\d+}/test.parquet",
#    run:
#        k = int(wildcards.k)
#        V = pd.read_parquet(input[0])
#        # converted to string to avoid error while matching
#        V["super_proximal"] = (V.tss_dist < 100).astype(str)
#        cols = ["super_proximal"]
#        V = match_cols(V[V.label], V[~V.label], cols, k=k, minimize_dist_col="tss_dist")
#        V = V.drop(columns=["super_proximal"])
#        print(V)
#        print(V.label.sum())
#        print(V.label.mean(), average_precision_score(V.label, -V.tss_dist))
#        V.to_parquet(output[0], index=False)
#

rule mendelian_dataset_v20:
    input:
        "results/omim/variants.annot_with_cre.annot_MAF.parquet",
        "results/gnomad/common.parquet",
        "results/tss.parquet",
    output:
        "results/dataset/mendelian_traits_v20_matched_{k,\d+}/test.parquet",
    run:
        k = int(wildcards.k)
        pos = pd.read_parquet(input[0])
        pos.maf = pos.maf.fillna(0)
        pos = pos[pos.maf < 0.1 / 100]
        pos = pos.drop(columns=["maf"])
        pos["label"] = True
        neg = pd.read_parquet(input[1])
        neg["label"] = False
        V = pd.concat([pos, neg], ignore_index=True)
        V = V[V.consequence.isin(TARGET_CONSEQUENCES)]
        assert len(V) == len(V.drop_duplicates(COORDINATES))

        V["start"] = V.pos - 1
        V["end"] = V.pos
        tss = pd.read_parquet(input[2], columns=["chrom", "start", "end"])
        V = bf.closest(V, tss).rename(columns={
            "distance": "tss_dist"
        }).drop(columns=["start", "end", "chrom_", "start_", "end_"])

        match_features = ["tss_dist"]

        consequences = V[V.label].consequence.unique()
        V_cs = []
        for c in consequences:
            print(c)
            V_c = V[V.consequence == c].copy()
            for f in match_features:
                V_c[f"{f}_scaled"] = RobustScaler().fit_transform(V_c[f].values.reshape(-1, 1))
            print(V_c.label.value_counts())
            V_c = match_columns_k(V_c, "label", [f"{f}_scaled" for f in match_features], k)
            V_c["match_group"] = c + "_" + V_c.match_group.astype(str)
            print(V_c.label.value_counts())
            print(V_c.groupby("label")[match_features].median())
            V_c.drop(columns=[f"{f}_scaled" for f in match_features], inplace=True)
            V_cs.append(V_c)
        V = pd.concat(V_cs, ignore_index=True)

        V = sort_variants(V)
        print(V)
        V.to_parquet(output[0], index=False)
    

rule mendelian_dataset_v21:
    input:
        "results/omim/variants.annot_with_cre.annot_MAF.parquet",
        "results/gnomad/common.parquet",
        "results/tss.parquet",
    output:
        "results/dataset/mendelian_traits_v21_matched_{k,\d+}/test.parquet",
    run:
        k = int(wildcards.k)
        pos = pd.read_parquet(input[0])
        pos.maf = pos.maf.fillna(0)
        pos = pos[pos.maf < 1 / 100]
        pos = pos.drop(columns=["maf"])
        pos["label"] = True
        neg = pd.read_parquet(input[1])
        neg["label"] = False
        V = pd.concat([pos, neg], ignore_index=True)
        V = V[V.consequence.isin(TARGET_CONSEQUENCES)]
        assert len(V) == len(V.drop_duplicates(COORDINATES))

        V["start"] = V.pos - 1
        V["end"] = V.pos
        tss = pd.read_parquet(input[2], columns=["chrom", "start", "end"])
        V = bf.closest(V, tss).rename(columns={
            "distance": "tss_dist"
        }).drop(columns=["start", "end", "chrom_", "start_", "end_"])

        match_features = ["tss_dist"]

        consequences = V[V.label].consequence.unique()
        V_cs = []
        for c in consequences:
            print(c)
            V_c = V[V.consequence == c].copy()
            for f in match_features:
                V_c[f"{f}_scaled"] = RobustScaler().fit_transform(V_c[f].values.reshape(-1, 1))
            print(V_c.label.value_counts())
            V_c = match_columns_k(V_c, "label", [f"{f}_scaled" for f in match_features], k)
            V_c["match_group"] = c + "_" + V_c.match_group.astype(str)
            print(V_c.label.value_counts())
            print(V_c.groupby("label")[match_features].median())
            V_c.drop(columns=[f"{f}_scaled" for f in match_features], inplace=True)
            V_cs.append(V_c)
        V = pd.concat(V_cs, ignore_index=True)

        V = sort_variants(V)
        print(V)
        V.to_parquet(output[0], index=False)