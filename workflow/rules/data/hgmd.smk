# version from Sei paper
rule hgmd_download:
    output:
        "results/hgmd/pathogenic_mutations/sei_data/hgmd.raw.matrix.csv",
    shell:
        """
        cd results/hgmd &&
        wget https://sei-files.s3.amazonaws.com/pathogenic_mutations.tar.gz &&
        tar -xzvf pathogenic_mutations.tar.gz
        """


rule hgmd_process:
    input:
        "results/hgmd/pathogenic_mutations/sei_data/hgmd.raw.matrix.csv",
        "results/genome.fa.gz",
    output:
        "results/hgmd/variants.parquet",
    run:
        V = pd.read_csv(
            input[0], usecols=["CHROM_hg19", "POS_hg19", "REF", "ALT", "phenotype"],
        )
        V = V.rename(columns={
            "CHROM_hg19": "chrom", "POS_hg19": "pos", "REF": "ref", "ALT": "alt",
            "phenotype": "trait",
        })
        print(V)
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
