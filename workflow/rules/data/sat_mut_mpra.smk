sat_mut_mpra_elements = [
    "F9",
    "GP1BA",
    "HBB",
    "HBG1",
    "HNF4A",
    "IRF4",
    "IRF6",
    "LDLR",
    "MSMB",
    "MYCrs6983267",
    "PKLR",
    "SORT1",
    "TCF7L2",
    "TERT",
    "ZFAND3",
    "ZRSh13"
]


rule sat_mut_mpra_download:
    output:
        temp("results/sat_mut_mpra/{element}/variants.vcf.gz"),
    wildcard_constraints:
        element="|".join(sat_mut_mpra_elements),
    shell:
        "wget -O {output} https://kircherlab.bihealth.org/download/CADD-development/v1.7/validation/regseq/SatMut.all.{wildcards.element}.vcf.gz"


rule sat_mut_mpra_process:
    input:
        "results/sat_mut_mpra/{element}/variants.vcf.gz",
    output:
        "results/sat_mut_mpra/{element}/test.parquet",
    wildcard_constraints:
        element="|".join(sat_mut_mpra_elements),
    run:
        V = pd.read_csv(
            input[0], sep="\t", dtype={"chrom": "str"}, comment="#", header=None,
            names=["chrom", "pos", "id", "ref", "alt", "qual", "filter", "INFO"],
            usecols=["chrom", "pos", "ref", "alt", "INFO"],
        )
        V["effect_size"] = V.INFO.str.extract(r"EF=([^;]+)").astype(float)
        V["p_value"] = V.INFO.str.extract(r"PV=([^;]+)").astype(float)
        V["barcodes"] = V.INFO.str.extract(r"BC=([^;]+)").astype(int)
        V.drop(columns=["INFO"], inplace=True)
        V = V.query("barcodes >= 10")
        V = V[V.ref.isin(NUCLEOTIDES) & V.alt.isin(NUCLEOTIDES)]
        V["label"] = V.effect_size.abs()
        V.to_parquet(output[0], index=False)


rule sat_mut_mpra_readme:
    input:
        expand("results/sat_mut_mpra/{element}/test.parquet", element=sat_mut_mpra_elements),
    output:
        "results/sat_mut_mpra/README.md",
    run:
        configs = []
        for e in sat_mut_mpra_elements:
            configs.append(dict(
                config_name=e,
                data_files=[dict(split="test", path=f"{e}/test.parquet")],
            ))
        metadata = dict(configs=configs)

        with open(output[0], "w") as readme:
            readme.write("---\n")
            readme.write(yaml.dump(metadata))
            readme.write("---\n")
