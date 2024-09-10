rule dataset_subset_all:
    input:
        "results/dataset/{dataset}/test.parquet",
    output:
        "results/dataset/{dataset}/subset/all.parquet",
    run:
        V = pd.read_parquet(input[0])
        V[COORDINATES].to_parquet(output[0], index=False)


rule dataset_subset_nonexonic:
    input:
        "results/dataset/{dataset}/test.parquet",
    output:
        "results/dataset/{dataset}/subset/nonexonic.parquet",
    run:
        V = pd.read_parquet(input[0])
        V = V[V.consequence.isin(NON_EXONIC + cre_classes + cre_flank_classes)]
        V[COORDINATES].to_parquet(output[0], index=False)


rule dataset_subset_consequence:
    input:
        "results/dataset/{dataset}/test.parquet",
    output:
        "results/dataset/{dataset}/subset/{c}.parquet",
    wildcard_constraints:
        c="|".join(other_consequences),
    run:
        V = pd.read_parquet(input[0])
        V = V[V.consequence == wildcards.c]
        V[COORDINATES].to_parquet(output[0], index=False)


rule dataset_subset_trait:
    input:
        "results/dataset/{dataset}/test.parquet",
    output:
        "results/dataset/{dataset}/subset/{trait}.parquet",
    wildcard_constraints:
        trait="|".join(traits_high_n),
    run:
        V = pd.read_parquet(input[0])
        V.trait = V.trait.str.split(",")
        target_size = len(V[V.match_group==V.match_group.iloc[0]])
        V = V[(~V.label) | (V.trait.apply(lambda x: wildcards.trait in x))]
        match_group_size = V.match_group.value_counts() 
        match_groups = match_group_size[match_group_size == target_size].index
        V = V[V.match_group.isin(match_groups)]
        V[COORDINATES].to_parquet(output[0], index=False)


rule dataset_subset_traits_n50:
    input:
        "results/dataset/{dataset}/test.parquet",
    output:
        "results/dataset/{dataset}/subset/traits_n50.parquet",
    wildcard_constraints:
        trait="|".join(traits_high_n),
    run:
        traits_n50 = set([
            "Height",
            "MCV",
            "MCH",
            "Mono",
            "Plt",
            "eBMD",
            "HbA1c",
            "RBC",
            "ALP",
            "IGF1",
            "HDLC",
            "Eosino",
            "LDLC",
            "SHBG",
            "AG",
            "Lym",
            "Hb",
            "GGT",
            "eGFRcys",
            "ApoA",
            "WBC",
            "eGFR",
            "TP",
        ])
        V = pd.read_parquet(input[0])
        V.trait = V.trait.str.split(",").apply(set)
        target_size = len(V[V.match_group==V.match_group.iloc[0]])
        V = V[(~V.label) | (V.trait.apply(lambda x: len(x.intersection(traits_n50)) > 0))]
        match_group_size = V.match_group.value_counts() 
        match_groups = match_group_size[match_group_size == target_size].index
        V = V[V.match_group.isin(match_groups)]
        print(V)
        V[COORDINATES].to_parquet(output[0], index=False)


rule dataset_subset_tissue:
    input:
        "results/dataset/{dataset}/test.parquet",
    output:
        "results/dataset/{dataset}/subset/{tissue}.parquet",
    wildcard_constraints:
        tissue="|".join(tissues),
    run:
        V = pd.read_parquet(input[0])
        V.tissue = V.tissue.str.split(",")
        target_size = len(V[V.match_group==V.match_group.iloc[0]])
        V = V[(~V.label) | (V.tissue.apply(lambda x: wildcards.tissue in x))]
        match_group_size = V.match_group.value_counts() 
        match_groups = match_group_size[match_group_size == target_size].index
        V = V[V.match_group.isin(match_groups)]
        V[COORDINATES].to_parquet(output[0], index=False)


rule run_classifier:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/dataset/{dataset}/subset/{subset}.parquet",
        lambda wildcards: expand("results/dataset/{{dataset}}/features/{features}.parquet", features=config["feature_sets"][wildcards.feature_set]),
    output:
        "results/dataset/{dataset}/preds/{subset}/{feature_set}.{classifier,LogisticRegression|RidgeRegression|RandomForest|XGBoost|PCALogisticRegression|FeatureSelectionLogisticRegression|BestFeature}.{split_mode,chrom|odd_even|pos_quarter}.parquet",
    threads:
        lambda wildcards: 1 if wildcards.classifier == "BestFeature" else workflow.cores
    run:
        V = pd.read_parquet(input[0])
        subset = pd.read_parquet(input[1])
        all_features = []
        for features, path in zip(config["feature_sets"][wildcards.feature_set], input[2:]):
            df = pd.read_parquet(path)
            df.columns = [f"{features}_{col}" for col in df.columns]
            all_features += df.columns.tolist()
            V = pd.concat([V, df], axis=1)
        V = subset.merge(V, on=COORDINATES, how="left")
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
        elif wildcards.split_mode == "pos_quarter":
            pos_subsets = np.array_split(np.unique(V.pos), 4)
            for pos_subset in pos_subsets:
                mask_train = V.pos.isin(pos_subset)
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
        "results/dataset/{dataset}/subset/{subset}.parquet",
        "results/dataset/{dataset}/features/{features}.parquet",
    output:
        "results/dataset/{dataset}/unsupervised_metrics/{subset}/{features}.csv",
    run:
        V = pd.read_parquet(input[0])
        subset = pd.read_parquet(input[1])
        df = pd.read_parquet(input[2])
        features = df.columns
        df = df.fillna(df.mean())
        V = pd.concat([V, df], axis=1)
        V = subset.merge(V, on=COORDINATES, how="left")
        balanced = V.label.sum() == len(V) // 2
        metric = roc_auc_score if balanced else average_precision_score
        metric_name = "AUROC" if balanced else "AUPRC"
        res = []
        for f in tqdm(features):
            score = max(metric(V.label, V[f]), metric(V.label, -V[f]))
            res.append([score, f])
        res = pd.DataFrame(res, columns=[metric_name, "feature"])
        res = res.sort_values(metric_name, ascending=False)
        res.to_csv(output[0], index=False)


#rule get_metrics:
#    input:
#        "results/dataset/{dataset}/test.parquet",
#        "results/dataset/{dataset}/preds/{model}.parquet",
#    output:
#        "results/dataset/{dataset}/metrics/{model}.csv",
#    run:
#        V = pd.read_parquet(input[0])
#        V["score"] = pd.read_parquet(input[1])["score"]
#        balanced = V.label.sum() == len(V) // 2
#        metric = roc_auc_score if balanced else average_precision_score
#        metric_name = "AUROC" if balanced else "AUPRC"
#        res = pd.DataFrame({
#            "Model": [wildcards.model],
#            metric_name: [metric(V.label, V.score)]
#        })
#        print(res)
#        res.to_csv(output[0], index=False)
#
#
#rule merge_metrics:
#    input:
#        lambda wildcards: expand(
#            "results/dataset/{{dataset}}/metrics/{model}.csv",
#            model=config["dataset_eval_models"][wildcards.dataset]
#        )
#    output:
#        "results/dataset/{dataset}/merged_metrics.csv",
#    run:
#        dfs = [pd.read_csv(f) for f in input]
#        df = pd.concat(dfs, axis=0)
#        df = df.sort_values(df.columns[-1], ascending=False)
#        print(df)
#        df.to_csv(output[0], index=False)
#
#
#rule plot_metrics:
#    input:
#        "results/dataset/{dataset}/test.parquet",
#        "results/dataset/{dataset}/merged_metrics.csv",
#    output:
#        "results/dataset/{dataset}/plot.pdf",
#    run:
#        V = pd.read_parquet(input[0])
#        n_pos, n_neg = V.label.sum(), len(V) - V.label.sum()
#        res = pd.read_csv(input[1])
#        metric = res.columns[-1]
#        if metric == "AUROC":
#            baseline = 0.5
#        elif metric == "AUPRC":
#            baseline = n_pos / (n_pos + n_neg)
#        plt.figure(figsize=(2,2))
#        g = sns.barplot(
#            data=res,
#            y="Model",
#            x=metric,
#            color="C0",
#        )
#        sns.despine()
#        sample_size = f"n={format_number(n_pos)} vs. {format_number(n_neg)}"
#        g.set(xlim=baseline, ylabel="")
#        title = f"{wildcards.dataset.split('/')[-1]}\n{sample_size}"
#        plt.title(title)
#        for bar, model in zip(g.patches, res.Model):
#            text = f'{bar.get_width():.3f}'
#
#            g.text(
#                max(bar.get_width(), baseline),  # X position, here at the end of the bar
#                bar.get_y() + bar.get_height()/2,  # Y position, in the middle of the bar
#                text,  # Text to be displayed, formatted to 3 decimal places
#                va='center'  # Vertical alignment
#            )
#        plt.savefig(output[0], bbox_inches="tight")


rule get_metrics_by_chrom:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/dataset/{dataset}/subset/{subset}.parquet",
        "results/dataset/{dataset}/preds/{subset}/{model}.parquet",
    output:
        "results/dataset/{dataset}/metrics_by_chrom/{subset}/{model}.csv",
    run:
        V = pd.read_parquet(input[0])
        subset = pd.read_parquet(input[1])
        V = subset.merge(V, on=COORDINATES, how="left")
        V["score"] = pd.read_parquet(input[2])["score"]
        balanced = V.label.sum() == len(V) // 2
        metric = roc_auc_score if balanced else average_precision_score
        metric_name = "AUROC" if balanced else "AUPRC"
        res = []
        for chrom in V.chrom.unique():
            V_chrom = V[V.chrom == chrom]
            res.append([chrom, len(V_chrom), wildcards.model, metric(V_chrom.label, V_chrom.score)])
        res = pd.DataFrame(res, columns=["chrom", "n", "Model", metric_name])
        res.to_csv(output[0], index=False)


rule get_metrics_by_chrom_weighted_average:
    input:
        "results/dataset/{dataset}/metrics_by_chrom/{subset}/{model}.csv",
    output:
        "results/dataset/{dataset}/metrics_by_chrom_weighted_average/{subset}/{model}.csv",
    run:
        res = pd.read_csv(input[0], dtype={"chrom": str})
        metric = res.columns[-1]
        res["weight"] = res.n / res.n.sum()
        res = pd.DataFrame({
            "Model": [wildcards.model],
            metric: [(res[metric] * res.weight).sum()]
        })
        res.to_csv(output[0], index=False)


#rule merge_metrics_by_chrom:
#    input:
#        lambda wildcards: expand(
#            "results/dataset/{{dataset}}/metrics_by_chrom/{model}.csv",
#            model=config["dataset_eval_models"][wildcards.dataset]
#        )
#    output:
#        "results/dataset/{dataset}/merged_metrics_by_chrom.csv",
#    run:
#        dfs = [pd.read_csv(f) for f in input]
#        df = pd.concat(dfs, axis=0)
#        print(df)
#        df.to_csv(output[0], index=False)
#
#
#rule plot_metrics_by_chrom:
#    input:
#        "results/dataset/{dataset}/test.parquet",
#        "results/dataset/{dataset}/merged_metrics_by_chrom.csv",
#    output:
#        "results/dataset/{dataset}/plot_by_chrom.pdf",
#    run:
#        V = pd.read_parquet(input[0])
#        n_pos, n_neg = V.label.sum(), len(V) - V.label.sum()
#        res = pd.read_csv(input[1])
#        metric = res.columns[-1]
#        if metric == "AUROC":
#            baseline = 0.5
#        elif metric == "AUPRC":
#            baseline = n_pos / (n_pos + n_neg)
#        plt.figure(figsize=(2,2))
#        g = sns.boxplot(
#            data=res,
#            y="Model",
#            x=metric,
#            color="C0",
#            order=res.groupby("Model")[metric].median().sort_values(ascending=False).index,
#        )
#        sns.despine()
#        sample_size = f"n={format_number(n_pos)} vs. {format_number(n_neg)}"
#        g.set(
#            #xlim=baseline,
#            ylabel=""
#        )
#        title = f"{wildcards.dataset.split('/')[-1]}\n{sample_size}"
#        plt.title(title)
#        plt.savefig(output[0], bbox_inches="tight")
#


rule get_metrics_by_odd_even:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/dataset/{dataset}/subset/{subset}.parquet",
        "results/dataset/{dataset}/preds/{subset}/{model}.parquet",
    output:
        "results/dataset/{dataset}/metrics_by_odd_even/{subset}/{model}.csv",
    run:
        V = pd.read_parquet(input[0])
        subset = pd.read_parquet(input[1])
        V = subset.merge(V, on=COORDINATES, how="left")
        V["score"] = pd.read_parquet(input[2])["score"]
        balanced = V.label.sum() == len(V) // 2
        metric = roc_auc_score if balanced else average_precision_score
        metric_name = "AUROC" if balanced else "AUPRC"
        res = []
        for name, chroms in zip(["odd", "even"], ODD_EVEN_CHROMS):
            V_chroms = V[V.chrom.isin(chroms)]
            res.append([name, len(V_chroms), wildcards.model, metric(V_chroms.label, V_chroms.score)])
        res = pd.DataFrame(res, columns=["chroms", "n", "Model", metric_name])
        res.to_csv(output[0], index=False)


rule get_metrics_by_odd_even_weighted_average:
    input:
        "results/dataset/{dataset}/metrics_by_odd_even/{subset}/{model}.csv",
    output:
        "results/dataset/{dataset}/metrics_by_odd_even_weighted_average/{subset}/{model}.csv",
    run:
        res = pd.read_csv(input[0])
        metric = res.columns[-1]
        res["weight"] = res.n / res.n.sum()
        res = pd.DataFrame({
            "Model": [wildcards.model],
            metric: [(res[metric] * res.weight).sum()]
        })
        res.to_csv(output[0], index=False)


rule get_metrics_by_pos_quarter:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/dataset/{dataset}/subset/{subset}.parquet",
        "results/dataset/{dataset}/preds/{subset}/{model}.parquet",
    output:
        "results/dataset/{dataset}/metrics_by_pos_quarter/{subset}/{model}.csv",
    run:
        V = pd.read_parquet(input[0])
        subset = pd.read_parquet(input[1])
        V = subset.merge(V, on=COORDINATES, how="left")
        V["score"] = pd.read_parquet(input[2])["score"]
        metric = lambda y_true, y_pred: pearsonr(y_true, y_pred)[0]
        metric_name = "Pearson"
        pos_subsets = np.array_split(np.unique(V.pos), 4)
        res = []
        for name, pos_subset in zip(range(len(pos_subsets)), pos_subsets):
            V_pos = V[V.pos.isin(pos_subset)]
            res.append([name, len(V_pos), wildcards.model, metric(V_pos.label, V_pos.score)])
        res = pd.DataFrame(res, columns=["pos_quarter", "n", "Model", metric_name])
        res.to_csv(output[0], index=False)


rule get_metrics_by_pos_quarter_weighted_average:
    input:
        "results/dataset/{dataset}/metrics_by_pos_quarter/{subset}/{model}.csv",
    output:
        "results/dataset/{dataset}/metrics_by_pos_quarter_weighted_average/{subset}/{model}.csv",
    run:
        res = pd.read_csv(input[0])
        metric = res.columns[-1]
        res["weight"] = res.n / res.n.sum()
        res = pd.DataFrame({
            "Model": [wildcards.model],
            metric: [(res[metric] * res.weight).sum()]
        })
        res.to_csv(output[0], index=False)