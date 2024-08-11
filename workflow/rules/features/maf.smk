rule download_gnomad_db:
    output:
        directory("results/gnomad_db"),
    run:
        from gnomad_db.database import gnomAD_DB

        url = "https://zenodo.org/record/6818606/files/gnomad_db_v3.1.2.sqlite3.gz?download=1"
        os.makedirs(output[0], exist_ok=True)
        gnomAD_DB.download_and_unzip(url, output[0])


rule maf_features:
    input:
        "results/gnomad_db",
    output:
        "results/features/{dataset}/minus_M{col,AF|AF_popmax}.parquet",
    threads:
        workflow.cores  # there seems to be an issue with concurrent db access
    run:
        from gnomad_db.database import gnomAD_DB

        V = load_dataset(wildcards["dataset"], split="test").to_pandas()
        db = gnomAD_DB(input[0], gnomad_version="v3")
        V["score"] = db.get_info_from_df(V, wildcards.col)[wildcards.col]
        V.score = V.score.where(V.score < 0.5, 1 - V.score)
        V.score = -V.score
        print(V.score.isna().sum())
        V = V.fillna(V["score"].mean())
        V[["score"]].to_parquet(output[0], index=False)
