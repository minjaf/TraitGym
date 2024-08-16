rule run_classifier:
    input:
        "results/dataset/{dataset}/test.parquet",
        lambda wildcards: expand("results/dataset/{{dataset}}/features/{features}.parquet", features=config["feature_sets"][wildcards.feature_set]),
    output:
        "results/dataset/{dataset}/preds/{feature_set}.{classifier,LogisticRegression|RandomForest|XGBoost}.{split_mode,chrom|odd_even}.parquet",
    threads:
        workflow.cores
    run:
        V = pd.read_parquet(input[0])
        all_features = []
        for features, path in zip(config["feature_sets"][wildcards.feature_set], input[1:]):
            df = pd.read_parquet(path)
            df.columns = [f"{features}_{col}" for col in df.columns]
            all_features += df.columns.tolist()
            V = pd.concat([V, df], axis=1)
        print(V)

        mask_train_list = []

        if wildcards.split_mode == "chrom":
            for chrom in V.chrom.unique():
                mask_train = V.chrom != chrom
                mask_train_list.append(mask_train)
        elif wildcards.split_mode == "odd_even":
            for chroms in ODD_EVEN_CHROMS:
                mask_train = V.chrom.isin(chroms)
                mask_train_list.append(mask_train)

        for mask_train in tqdm(mask_train_list):
            mask_test = ~mask_train
            V.loc[mask_test, "score"] = classifier_map[wildcards.classifier](
                V[mask_train], V[mask_test], all_features
            )

        V[["score"]].to_parquet(output[0], index=False)


rule eval_unsupervised_features:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/dataset/{dataset}/features/{features}.parquet",
    output:
        "results/dataset/{dataset}/unsupervised_metrics/{features}.csv",
    run:
        V = pd.read_parquet(input[0])
        features = pd.read_parquet(input[1])
        features = features.fillna(features.mean())
        balanced = V.label.sum() == len(V) // 2
        metric = roc_auc_score if balanced else average_precision_score
        metric_name = "AUROC" if balanced else "AUPRC"
        res = []
        for f in tqdm(features.columns):
            score = max(metric(V.label, features[f]), metric(V.label, -features[f]))
            res.append([score, f])
        res = pd.DataFrame(res, columns=[metric_name, "feature"])
        res = res.sort_values(metric_name, ascending=False)
        res.to_csv(output[0], index=False)


#rule unsupervised_l2_score:
#    output:
#        "results/preds/{dataset}/{features}.Unsupervised.L2.parquet",
#    run:
#        df = pd.read_parquet(f"https://huggingface.co/datasets/{wildcards.dataset}/resolve/main/features/{wildcards.features}.parquet")
#        df = df.fillna(df.mean())
#        df["score"] = np.linalg.norm(df, axis=1)
#        df[["score"]].to_parquet(output[0], index=False)
#
#
#rule unsupervised_scalar_score:
#    output:
#        "results/preds/{dataset}/{features}.Unsupervised.scalar.parquet",
#    run:
#        df = pd.read_parquet(f"https://huggingface.co/datasets/{wildcards.dataset}/resolve/main/features/{wildcards.features}.parquet")
#        assert df.shape[1] == 1
#        df = df.fillna(df.mean())
#        df = df.rename(columns={df.columns[0]: "score"})
#        df.to_parquet(output[0], index=False)


rule get_metrics:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/dataset/{dataset}/preds/{model}.parquet",
    output:
        "results/dataset/{dataset}/metrics/{model}.csv",
    run:
        V = pd.read_parquet(input[0])
        V["score"] = pd.read_parquet(input[1])["score"]
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
            "results/dataset/{{dataset}}/metrics/{model}.csv",
            model=config["dataset_eval_models"][wildcards.dataset]
        )
    output:
        "results/dataset/{dataset}/merged_metrics.csv",
    run:
        dfs = [pd.read_csv(f) for f in input]
        df = pd.concat(dfs, axis=0)
        df = df.sort_values(df.columns[-1], ascending=False)
        print(df)
        df.to_csv(output[0], index=False)


rule plot_metrics:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/dataset/{dataset}/merged_metrics.csv",
    output:
        "results/dataset/{dataset}/plot.pdf",
    run:
        V = pd.read_parquet(input[0])
        n_pos, n_neg = V.label.sum(), len(V) - V.label.sum()
        res = pd.read_csv(input[1])
        metric = res.columns[-1]
        if metric == "AUROC":
            baseline = 0.5
        elif metric == "AUPRC":
            baseline = n_pos / (n_pos + n_neg)
        plt.figure(figsize=(2,2))
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


rule get_metrics_by_chrom:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/dataset/{dataset}/preds/{model}.parquet",
    output:
        "results/dataset/{dataset}/metrics_by_chrom/{model}.csv",
    run:
        V = pd.read_parquet(input[0])
        V["score"] = pd.read_parquet(input[1])["score"]
        balanced = V.label.sum() == len(V) // 2
        metric = roc_auc_score if balanced else average_precision_score
        metric_name = "AUROC" if balanced else "AUPRC"
        res = []
        for chrom in V.chrom.unique():
            V_chrom = V[V.chrom == chrom]
            res.append([chrom, wildcards.model, metric(V_chrom.label, V_chrom.score)])
        res = pd.DataFrame(res, columns=["chrom", "Model", metric_name])
        print(res)
        res.to_csv(output[0], index=False)


rule merge_metrics_by_chrom:
    input:
        lambda wildcards: expand(
            "results/dataset/{{dataset}}/metrics_by_chrom/{model}.csv",
            model=config["dataset_eval_models"][wildcards.dataset]
        )
    output:
        "results/dataset/{dataset}/merged_metrics_by_chrom.csv",
    run:
        dfs = [pd.read_csv(f) for f in input]
        df = pd.concat(dfs, axis=0)
        print(df)
        df.to_csv(output[0], index=False)


rule plot_metrics_by_chrom:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/dataset/{dataset}/merged_metrics_by_chrom.csv",
    output:
        "results/dataset/{dataset}/plot_by_chrom.pdf",
    run:
        V = pd.read_parquet(input[0])
        n_pos, n_neg = V.label.sum(), len(V) - V.label.sum()
        res = pd.read_csv(input[1])
        metric = res.columns[-1]
        if metric == "AUROC":
            baseline = 0.5
        elif metric == "AUPRC":
            baseline = n_pos / (n_pos + n_neg)
        plt.figure(figsize=(2,2))
        g = sns.boxplot(
            data=res,
            y="Model",
            x=metric,
            color="C0",
            order=res.groupby("Model")[metric].median().sort_values(ascending=False).index,
        )
        sns.despine()
        sample_size = f"n={format_number(n_pos)} vs. {format_number(n_neg)}"
        g.set(
            #xlim=baseline,
            ylabel=""
        )
        title = f"{wildcards.dataset.split('/')[-1]}\n{sample_size}"
        plt.title(title)
        plt.savefig(output[0], bbox_inches="tight")
