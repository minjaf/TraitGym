insights_metadata = pl.read_excel(
    'config/insights.xlsx', sheet_name='ST2_overview_traits',
    read_options=dict(header_row=2)
).to_pandas()
finemap_methods = ["susie", "finemap"]


rule insights_finngen_download:
    output:
        "results/insights/finngen/raw/{method}/{trait}.parquet",
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
                columns=["chromosome", "position", "allele1", "allele2", "maf", "p", "prob"]
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

    
rule insights_finngen_avg_methods:
    input:
        expand(
            "results/insights/finngen/raw/{method}/{{trait}}.parquet",
            method=finemap_methods
        ),
    output:
        "results/insights/finngen/raw/avg/{trait}.parquet",
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
