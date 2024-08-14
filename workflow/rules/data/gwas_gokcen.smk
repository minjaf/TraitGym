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


rule gwas_gokcen_experiment1:
    input:
        "results/gwas_gokcen/processed/{trait}.parquet",
        expand("results/intervals/cre_{c}.parquet", c=cre_classes),
    output:
        "results/gwas_gokcen/experiment1/{trait}.parquet",
    run:
        pass
