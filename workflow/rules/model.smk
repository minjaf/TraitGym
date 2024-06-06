rule benchmark_set:
    output:
        "results/benchmark_set/{dataset}.parquet",
    run:
        V = load_dataset(wildcards.dataset, split="test").to_pandas()
        print(V)
        all_features = []
        for features in config["baseline_features"][wildcards.dataset]:
            df = pd.read_parquet(f"https://huggingface.co/datasets/{wildcards.dataset}/resolve/main/features/{features}.parquet")
            df.columns = [f"{features}_{col}" for col in df.columns]
            all_features += df.columns.tolist()
            V = pd.concat([V, df], axis=1)
        V = V.dropna(subset=all_features)
        if "match_group" in V.columns:
            V = V[V.duplicated("match_group", keep=False)]
        print(V)
        V[COORDINATES].to_parquet(output[0], index=False)


rule run_classifier:
    input:
        "results/benchmark_set/{dataset}.parquet",
    output:
        "results/preds/{dataset}/{feature_set}.{classifier,LogisticRegression|XGBoost}.parquet",
    threads:
        workflow.cores
    run:
        V_full = load_dataset(wildcards.dataset, split="test").to_pandas()
        print(V_full)
        all_features = []
        for features in config["feature_sets"][wildcards.feature_set]:
            df = pd.read_parquet(f"https://huggingface.co/datasets/{wildcards.dataset}/resolve/main/features/{features}.parquet")
            df.columns = [f"{features}_{col}" for col in df.columns]
            all_features += df.columns.tolist()
            V_full = pd.concat([V_full, df], axis=1)
        print(V_full)
        V_subset = pd.read_parquet(input[0])
        V = V_full.merge(V_subset, on=COORDINATES, how="inner")
        print(V)

        #for chrom in tqdm(V.chrom.unique()):
        #    mask_train = V.chrom != chrom
        for chroms in tqdm(ODD_EVEN_CHROMS):
            mask_train = V.chrom.isin(chroms)
            mask_test = ~mask_train
            V.loc[mask_test, "score"] = classifier_map[wildcards.classifier](
                V[mask_train], V[mask_test], all_features
            )

        V_full[COORDINATES].merge(
            V[COORDINATES + ["score"]], on=COORDINATES, how="left"
        ).to_parquet(output[0], index=False)


# let's just do L2 norm for now
rule unsupervised_score:
    output:
        "results/preds/{dataset}/{features}.Unsupervised.parquet",
    run:
        df = pd.read_parquet(f"https://huggingface.co/datasets/{wildcards.dataset}/resolve/main/features/{wildcards.features}.parquet")
        df["score"] = np.linalg.norm(df, axis=1)
        df[["score"]].to_parquet(output[0], index=False)


rule get_metrics:
    input:
        "results/preds/{dataset}/{model}.parquet",
        "results/benchmark_set/{dataset}.parquet",
    output:
        "results/metrics/{dataset}/{model}.csv",
    run:
        V = load_dataset(wildcards.dataset, split="test").to_pandas()
        V["score"] = pd.read_parquet(input[0])["score"]
        V = V.merge(pd.read_parquet(input[1]), on=COORDINATES, how="inner")
        balanced = V.label.sum() == len(V) // 2
        metric = roc_auc_score if balanced else average_precision_score
        metric_name = "AUROC" if balanced else "AUPRC"
        res = pd.DataFrame({
            "Model": [wildcards.model],
            metric_name: [metric(V.label, V.score)]
        })
        print(res)
        res.to_csv(output[0], index=False)


rule merge_metrics:
    input:
        lambda wildcards: expand(
            "results/metrics/{{dataset}}/{model}.csv",
            model=config["dataset_eval_models"][wildcards.dataset]
        )
    output:
        "results/merged_metrics/{dataset}.csv",
    run:
        dfs = [pd.read_csv(f) for f in input]
        df = pd.concat(dfs, axis=0)
        df = df.sort_values(df.columns[-1], ascending=False)
        print(df)
        df.to_csv(output[0], index=False)


rule plot_metrics:
    input:
        "results/benchmark_set/{dataset}.parquet",
        "results/merged_metrics/{dataset}.csv",
    output:
        "results/plots/{dataset}.pdf",
    run:
        V = load_dataset(wildcards.dataset, split="test").to_pandas()
        V = V.merge(pd.read_parquet(input[0]), on=COORDINATES, how="inner")
        n_pos, n_neg = V.label.sum(), len(V) - V.label.sum()
        res = pd.read_csv(input[1])
        metric = res.columns[-1]
        if metric == "AUROC":
            baseline = 0.5
        elif metric == "AUPRC":
            baseline = n_pos / (n_pos + n_neg)
        plt.figure(figsize=(3,3))
        g = sns.barplot(
            data=res,
            y="Model",
            x=metric,
            color="C0",
        )
        sns.despine()
        sample_size = f"n={format_number(n_pos)} vs. {format_number(n_neg)}"
        g.set(xlim=baseline, ylabel="")
        title = f"{wildcards.dataset.split('/')[-1]}\n{sample_size}"
        plt.title(title)
        for bar, model in zip(g.patches, res.Model):
            text = f'{bar.get_width():.3f}'

            g.text(
                max(bar.get_width(), baseline),  # X position, here at the end of the bar
                bar.get_y() + bar.get_height()/2,  # Y position, in the middle of the bar
                text,  # Text to be displayed, formatted to 3 decimal places
                va='center'  # Vertical alignment
            )
        plt.savefig(output[0], bbox_inches="tight")