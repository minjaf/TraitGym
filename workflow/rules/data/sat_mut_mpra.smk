sat_mut_mpra_promoters = [
    "F9",
    "GP1BA",
    "HBB",
    "HBG1",
    "HNF4A",
    "LDLR",
    "MSMB",
    "PKLR",
    "TERT",
]

sat_mut_mpra_enhancers = [
    "IRF4",
    "IRF6",
    "MYCrs6983267",
    "SORT1",
    "TCF7L2",
    "ZFAND3",
    "ZRSh13"
]

sat_mut_mpra_proximity_mapping = {
    "proximal": sat_mut_mpra_promoters,
    "distal": sat_mut_mpra_enhancers,
}

sat_mut_mpra_elements = sat_mut_mpra_promoters + sat_mut_mpra_enhancers


rule sat_mut_mpra_download:
    output:
        temp("results/sat_mut_mpra/{element}.vcf.gz"),
    wildcard_constraints:
        element="|".join(sat_mut_mpra_elements),
    shell:
        "wget -O {output} https://kircherlab.bihealth.org/download/CADD-development/v1.7/validation/regseq/SatMut.all.{wildcards.element}.vcf.gz"


rule sat_mut_mpra_process:
    input:
        "results/sat_mut_mpra/{element}.vcf.gz",
    output:
        "results/sat_mut_mpra/{element}.parquet",
    wildcard_constraints:
        element="|".join(sat_mut_mpra_elements),
    run:
        V = pd.read_csv(
            input[0], sep="\t", dtype={"chrom": "str"}, comment="#", header=None,
            names=["chrom", "pos", "id", "ref", "alt", "qual", "FILTER", "INFO"],
            usecols=["chrom", "pos", "ref", "alt", "FILTER", "INFO"],
        )
        V = V[V.FILTER=="SIGN"]
        V["label"] = V.INFO.str.extract(r"EF=([^;]+)").astype(float)
        V["element"] = wildcards.element
        V.drop(columns=["FILTER", "INFO"], inplace=True)
        V = V[V.ref.isin(NUCLEOTIDES) & V.alt.isin(NUCLEOTIDES)]
        V.to_parquet(output[0], index=False)


rule sat_mut_mpra_merge:
    input:
        lambda wildcards: expand("results/sat_mut_mpra/{element}.parquet", element=sat_mut_mpra_proximity_mapping[wildcards.proximity]),
    output:
        "results/dataset/sat_mut_mpra_{proximity}/test.parquet",
    run:
        V = pd.concat([pd.read_parquet(i) for i in input], ignore_index=True)
        V["effect_size"] = V.label
        V.label = V.label.abs()
        V = sort_chrom_pos(V)
        print(V)
        V.to_parquet(output[0], index=False)


#rule sat_mut_mpra_subset_element:
#    input:
#        "results/dataset/sat_mut_mpra/test.parquet",
#    output:
#        "results/dataset/sat_mut_mpra/subset/{e}.parquet",
#    wildcard_constraints:
#        e="|".join(sat_mut_mpra_elements),
#    run:
#        V = pd.read_parquet(input[0])
#        V = V[V.element == wildcards.e]
#        V[COORDINATES].to_parquet(output[0], index=False)
