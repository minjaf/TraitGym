evo2_versions = list(config["evo2"].keys())

n_shards = 6
shards = np.arange(n_shards)


#rule evo2_run_vep:
#    input:
#        "results/dataset/{dataset}/test.parquet",
#        "results/genome.fa.gz",
#    output:
#        "results/dataset/{dataset}/features/evo2_{version}.parquet",
#    wildcard_constraints:
#        version="|".join(evo2_versions),
#    params:
#        window_size=lambda wildcards: config["evo2"][wildcards.version]["window_size"],
#        model_path=lambda wildcards: config["evo2"][wildcards.version]["model_path"],
#        per_device_batch_size=lambda wildcards: config["evo2"][wildcards.version]["per_device_batch_size"],
#    threads:
#        8
#    priority: 999
#    shell:
#        """
#        python \
#        workflow/scripts/run_vep_evo2.py {input} {params.window_size} \
#        {params.model_path} {output} --is_file --dataloader_num_workers 8 \
#        --per_device_batch_size {params.per_device_batch_size} \
#        """


rule evo2_merge_shards:
    input:
        expand(
            "results/dataset/{{dataset}}/features/shard/{shard}/evo2_40b.parquet",
            shard=shards
        ),
    output:
        "results/dataset/{dataset}/features/evo2_40b.parquet",
    run:
        df = pd.concat([pd.read_parquet(path) for path in input])
        df.to_parquet(output[0], index=False)


rule evo2_extract_llr:
    input:
        "results/dataset/{dataset}/features/evo2_{version}.parquet",
    output:
        "results/dataset/{dataset}/features/evo2_{version}_LLR.parquet",
    wildcard_constraints:
        version="|".join(evo2_versions),
    run:
        (
            pl.read_parquet(input[0], columns=["llr"])
            .rename({"llr": "score"})
            .write_parquet(output[0])
        )


rule evo2_extract_embeddings:
    input:
        "results/dataset/{dataset}/features/evo2_{version}.parquet",
    output:
        "results/dataset/{dataset}/features/evo2_{version}_Embeddings.parquet",
    wildcard_constraints:
        version="|".join(evo2_versions),
    run:
        df = pl.read_parquet(input[0])
        print(df)
        df = df[:, 1:]
        print(df)
        df.write_parquet(output[0])
