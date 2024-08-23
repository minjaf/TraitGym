insights_metadata = pl.read_excel(
    'config/insights.xlsx', sheet_name='ST2_overview_traits',
    read_options=dict(header_row=2)
).to_pandas()
fg_traits = insights_metadata.query('cohort == "FG"').trait.values
fg_traits_exclude = [
    "Malignant_Neoplasms",
    "Urolithiasis",
    "Hypertension",
    "Glaucoma",
    "PAD",
    "RA",
    "Benign_Neoplasms",
    "Asthma",
    "Gout",
    "ICP",
    "UC",
]
fg_traits = [t for t in fg_traits if t not in fg_traits_exclude]
finemap_methods = ["susie", "finemap"]


# gbenegas: some traits are missing in the latest release
# so using my discretion to match traits
#trait_code_override = {
#    "Malignant_Neoplasms": "C3_CANCER_EXALLC",  # this excludes all cancers. I'm not sure if it's the same
#    #"Urolithiasis": pass,  # non-core endpoint, seems like it was subdivided, no finemapping
#    "Hypertension": "I9_HypTens",
#    "Glaucoma": "H7_GLAUCOMA",
#    #"PAD": None,  # non-core endpoint, seems like it was subdivided, no finemapping
#    #"RA": None,  # non-core endpoint, seems like it was subdivided, no finemapping
#    #"Benign_Neoplasms": None,  # seems like it was subdivided
#    #"Asthma": None,
#    # "Gout": None,
#    "ICP": None,
#    "UC": None,
#}


rule insights_fg_download:
    output:
        temp("results/insights/fg/raw/{method}/{trait}.parquet"),
    wildcard_constraints:
        method="|".join(finemap_methods),
    run:
        trait_code = (
            insights_metadata.query(f'cohort == "FG" and trait == "{wildcards.trait}"')
            .definition.str.split(" ").str[-1].iloc[0]
        )
        path = f"https://storage.googleapis.com/finngen-public-data-r11/finemap/full/{wildcards.method}/finngen_R11_{trait_code}.{wildcards.method.upper()}.snp.bgz"
        print(path)
        V = (
            pl.read_csv(
                path, separator='\t',
                columns=["chromosome", "position", "allele1", "allele2", "maf", "p", "prob"],
                schema_overrides={"prob": float},
            )
            .rename({
                "chromosome": "chrom",
                "position": "pos",
                "allele1": "ref",
                "allele2": "alt",
                "prob": "pip"
            })
            .select(["chrom", "pos", "ref", "alt", "maf", "p", "pip"])
        )
        print(V)
        V.write_parquet(output[0])

    
rule insights_fg_avg_methods:
    input:
        expand(
            "results/insights/fg/raw/{method}/{{trait}}.parquet",
            method=finemap_methods
        ),
    output:
        temp("results/insights/fg/raw/avg/{trait}.parquet"),
    run:
        V = pl.concat([pl.read_parquet(f) for f in input])
        V = (
            V.group_by(["chrom", "pos", "ref", "alt"])
            .agg(
                pl.first("maf"),
                pl.first("p"),
                pl.mean("pip"),
                (pl.max("pip") - pl.min("pip")).alias("pip_diff"),
            )
        )
        print(V)
        V = V.filter(pl.col("pip_diff") < 0.05).drop("pip_diff")
        print(V)
        V.write_parquet(output[0])


rule insights_fg_process:
    input:
        "results/insights/fg/raw/avg/{trait}.parquet",
        "results/genome.fa.gz",
    output:
        "results/insights/fg/processed/{trait}.parquet",
    run:
        V = pd.read_parquet(input[0])
        V.chrom = V.chrom.str.replace("chr", "")
        V = filter_chroms(V)
        V = filter_snp(V)
        V = sort_variants(V)
        genome = Genome(input[1])
        V = check_ref_alt(V, genome)
        V.to_parquet(output[0], index=False)


# this should merge all the 3 cohorts, in the future
rule insights_merge_coords:
    input:
        expand("results/insights/fg/processed/{trait}.parquet", trait=fg_traits),
    output:
        "results/insights/coords.parquet",
    run:
        V = (
            pl.concat([pl.read_parquet(f, columns=COORDINATES) for f in input])
            .unique()
            .to_pandas()
        )
        V = sort_variants(V)
        print(V)
        V.to_parquet(output[0], index=False)
