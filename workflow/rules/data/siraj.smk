rule siraj_download_st_other:
    output:
        "results/siraj/st_other.xlsx",
    shell:
        "wget -O {output} https://www.biorxiv.org/content/biorxiv/early/2024/05/06/2024.05.05.592437/DC9/embed/media-9.xlsx?download=true"


rule siraj_pip90:
    input:
        "results/siraj/st_other.xlsx",
        "results/genome.fa.gz",
    output:
        "results/siraj/{c,gwas|eqtl}_pip90/test.parquet",
    run:
        V = (
            pl.read_excel(
                input[0], sheet_name="S4.MPRA PrC", read_options=dict(header_row=1)
            )
            .to_pandas()
        )
        if wildcards.c == "gwas":
            V = V[V["Variant Type"].str.contains("Complex traits")]
        elif wildcards.c == "eqtl":
            V = V[V["Variant Type"].str.contains("eQTL")]
        V["label"] = V["Variant Type"].str.contains("test")
        V["chrom"] = V.Variant.str.split(":").str[0].str.replace("chr", "")
        V["pos"] = V.Variant.str.split(":").str[1].astype(int)
        V["ref"] = V.Variant.str.split(":").str[2]
        V["alt"] = V.Variant.str.split(":").str[3]
        V = V[["chrom", "pos", "ref", "alt", "label"]]
        V = V[(V.ref.isin(NUCLEOTIDES)) & (V.alt.isin(NUCLEOTIDES))]
        V = lift_hg19_to_hg38(V)
        V = V[V.pos != -1]
        genome = Genome(input[1])
        V = check_ref_alt(V, genome)
        V = (
            V.groupby("label")
            .apply(lambda x: x.sample(n=V.label.value_counts().min(), random_state=42))
            .reset_index(drop=True)
        )
        V = sort_chrom_pos(V)
        V.to_parquet(output[0], index=False)
