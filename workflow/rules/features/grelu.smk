rule grelu_features:
    output:
        "results/features/{dataset}/{model,Enformer|Borzoi}_L2.parquet",
    threads:
        workflow.cores
    run:
        import grelu.resources
        import grelu.variant

        V = load_dataset(wildcards.dataset, split="test").to_pandas()
        V.chrom = "chr" + V.chrom
        print(V)

        if wildcards.model == "Enformer":
            project, model_name = "enformer", "human"
        elif wildcards.model == "Borzoi":
            project, model_name = "borzoi", "human_fold0"
        model = grelu.resources.load_model(project=project, model_name=model_name)

        chrom_sizes = (
            bf.fetch_chromsizes("hg38").to_frame().reset_index()
            .rename(columns={"index": "chrom"})
        )
        L = model.data_params["train_seq_len"]
        V["start"] = V.pos - 1 - L // 2
        V["end"] = V.start + L
        V = pd.merge(V, chrom_sizes, on="chrom", how="left")
        V["within_bounds"] = (V.start >= 0) & (V.end <= V.length)
        print(f"{V.within_bounds.value_counts()=}")

        def log2fc_1p(x: np.ndarray, y: np.ndarray) -> np.ndarray:
            return np.log2(1+x) - np.log2(1+y)

        lfc = grelu.variant.predict_variant_effects(
            variants=V.loc[V.within_bounds, COORDINATES],
            model=model, 
            devices=list(range(torch.cuda.device_count())),
            num_workers=threads,
            batch_size=8,
            genome="hg38",
            compare_func=log2fc_1p,
            return_ad=False,
            # Reverse complement the ref/alt predictions and average them.
            rc=True,
        )
        l2_score = np.linalg.norm(lfc, axis=2)
        columns = model.data_params['tasks']["name"]
        V.loc[V.within_bounds, columns] = l2_score
        V[columns].to_parquet(output[0], index=False)
