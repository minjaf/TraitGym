rule dataset_to_vcf:
    output:
        "results/data/{dataset}.vcf.gz",
    run:
        V = load_dataset(wildcards.dataset, split="test").to_pandas()
        V["id"] = "."
        V[["chrom", "pos", "id", "ref", "alt"]].to_csv(
            output[0], sep="\t", index=False, header=False,
        )


# obtained from the website
rule cadd_process:
    input:
        "results/data/cadd/{dataset}.tsv.gz",
    output:
        "results/features/{dataset}/CADD_RawScore.parquet",
        "results/features/{dataset}/CADD_Annot.parquet",
    run:
        V = load_dataset(wildcards.dataset, split="test").to_pandas()
        df = pd.read_csv(input[0], sep="\t", skiprows=1, dtype={"#Chrom": "str"})
        df = df.rename(columns={"#Chrom": "chrom", "Pos": "pos", "Ref": "ref", "Alt": "alt"})
        
        # some variants have separate rows for AnnoType = Transcript vs. Intergenic
        # we want features that are independent of the AnnoType
        dup = df[df.duplicated(COORDINATES, keep=False)][COORDINATES].iloc[0]
        df2 = df[df[COORDINATES].eq(dup).all(axis=1)]
        independent_cols = [c for c in df2.columns if df2[c].nunique() == 1]
        df = df.drop_duplicates(COORDINATES)[independent_cols]

        df = df.loc[:, df.isna().mean() < 0.1]
        df = df.loc[:, df.nunique() > 1]
        # drop string columns
        for c in df.columns:
            if c in COORDINATES: continue
            if df[c].apply(lambda x: isinstance(x, str)).all():
                df.drop(columns=[c], inplace=True)
        # Fill NaN with the mean of the column
        other_cols = [c for c in df.columns if c not in COORDINATES]
        for c in other_cols:
            df[c] = df[c].fillna(df[c].mean())
        V = V.merge(df, how="left", on=COORDINATES)
        V[["RawScore"]].to_parquet(output[0], index=False)
        V[other_cols].to_parquet(output[1], index=False)
