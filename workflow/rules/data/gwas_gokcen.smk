gwas_gokcen_metadata = pd.read_csv(
    "results/gwas_gokcen/raw/disease_list.txt", delim_whitespace=True, header=None,
    names=["trait", "source", "name"]
)
gwas_gokcen_metadata["path"] = (
    gwas_gokcen_metadata.trait + "-" + gwas_gokcen_metadata.source
)
gwas_gokcen_metadata.set_index("trait", inplace=True)


def get_input_gwas_gokcen_process(wildcards):
    path = gwas_gokcen_metadata.loc[wildcards.trait, "path"]
    return f"results/gwas_gokcen/raw/{path}.susie.gwfinemap.b38.gz"


rule gwas_gokcen_process:
    input:
        get_input_gwas_gokcen_process,
        "results/genome.fa.gz",
    output:
        "results/gwas_gokcen/processed/{trait}.parquet",
    run:
        V = pl.read_csv(
            input[0], separator="\t",
            columns=["CHR", "BP", "A1", "A2", "MAF", "PIP"],
            new_columns=["chrom", "pos", "ref", "alt", "MAF", "PIP"],
        ).to_pandas()
        V.chrom = V.chrom.str.replace("chr", "")
        V = filter_chroms(V)
        V = filter_snp(V)
        genome = Genome(input[1])
        V = check_ref_alt(V, genome)
        V.to_parquet(output[0], index=False)


rule gwas_gokcen_experiment1_filt:
    input:
        "results/gwas_gokcen/processed/{trait}.parquet",
        expand("results/intervals/cre_{c}.parquet", c=cre_classes),
    output:
        "results/gwas_gokcen/experiment1/filt/{trait}.parquet",
    run:
        V = pd.read_parquet(input[0])
        print(V.shape)
        V["start"] = V.pos - 1
        V["end"] = V.pos
        for c, path in zip(cre_classes, input[1:]):
            I = pd.read_parquet(path)
            V = bf.coverage(V, I)
            V[c] = V.coverage > 0
            V = V.drop(columns=["coverage"])
        V = V.drop(columns=["start", "end"])
        V = V[V[cre_classes].any(axis=1)]
        V.loc[V.PIP > 0.9, "label"] = True
        V.loc[V.PIP < 0.01, "label"] = False
        V = V.dropna(subset=["label"])
        print(V)
        V.to_parquet(output[0], index=False)


rule gwas_gokcen_experiment1_match:
    input:
        "results/gwas_gokcen/experiment1/filt/{trait}.parquet",
        "results/tss.parquet",
    output:
        "results/dataset/gwas_gokcen_experiment1_{trait}/test.parquet",
    run:
        V = pd.read_parquet(input[0])

        V["start"] = V.pos - 1
        V["end"] = V.pos

        tss = pd.read_parquet(input[1], columns=["chrom", "start", "end"])

        V = bf.closest(V, tss).rename(columns={
            "distance": "tss_dist"
        }).drop(columns=["start", "end", "chrom_", "start_", "end_"])

        match_features = ["MAF", "tss_dist"]

        V_cs = []
        for c in cre_classes:
            print(c)
            V_c = V[V[c]].copy()
            for f in match_features:
                V_c[f"{f}_scaled"] = RobustScaler().fit_transform(V_c[f].values.reshape(-1, 1))
            print(V_c.label.value_counts())
            V_c = match_columns(V_c, "label", [f"{f}_scaled" for f in match_features])
            V_c["match_group"] = c + "_" + V_c.match_group.astype(str)
            print(V_c.label.value_counts())
            print(V_c.groupby("label")[match_features].median())
            V_c.drop(columns=[f"{f}_scaled" for f in match_features], inplace=True)
            V_cs.append(V_c)
        V = pd.concat(V_cs, ignore_index=True)
        V = sort_chrom_pos(V)
        print(V)
        V.to_parquet(output[0], index=False)
