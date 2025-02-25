n_shards = 8
shards = np.arange(n_shards)


rule dataset_shard:
    input:
        "results/dataset/{dataset}/test.parquet",
    output:
        expand("results/dataset/{{dataset}}/shard/{shard}.parquet", shard=shards),
    run:
        df = pd.read_parquet(input[0])
        for path, df_s in zip(output, np.array_split(df, n_shards)):
            df_s.to_parquet(path, index=False)


rule evo2_run_vep:
    input:
        "results/dataset/{dataset}/shard/{shard}.parquet",
        "results/genome.fa.gz",
    output:
        "results/dataset/{dataset}/features/shard/{shard}/evo2_{version}.parquet",
    resources:
        nvidia_gpu=1
    params:
        window_size=lambda wildcards: config["evo2"][wildcards.version]["window_size"],
        model_path=lambda wildcards: config["evo2"][wildcards.version]["model_path"],
        per_device_batch_size=lambda wildcards: config["evo2"][wildcards.version]["per_device_batch_size"],
    threads:
        8
    priority: 999
    shell:
        """
        python \
        workflow/scripts/run_vep_evo2.py {input} {params.window_size} \
        {params.model_path} {output} --is_file --dataloader_num_workers 8 \
        --per_device_batch_size {params.per_device_batch_size} \
        """

#        OMP_NUM_THREADS=8 \
#        torchrun --nproc_per_node $(echo $CUDA_VISIBLE_DEVICES | awk -F',' '{{print NF}}') \


rule evo2_merge_shards:
    input:
        expand(
            "results/dataset/{{dataset}}/features/shard/{shard}/evo2_{{version}}.parquet",
            shard=shards
        ),
    output:
        "results/dataset/{dataset}/features/evo2_{version}.parquet",
    run:
        df = pd.concat([pd.read_parquet(path) for path in input])
        df.to_parquet(output[0], index=False)
