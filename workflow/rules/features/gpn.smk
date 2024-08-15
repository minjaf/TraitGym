rule gpn_run_vep_embed_dists:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/genome.fa.gz",
    output:
        "results/dataset/{dataset}/features/GPN_EmbedDists.parquet",
    threads:
        workflow.cores
    shell:
        """
        python \
        -m gpn.ss.run_vep_embed_dists {input} {config[gpn][window_size]} \
        {config[gpn][model_path]} {output} --is_file \
        --per-device-batch-size 2048 --dataloader-num-workers {threads}
        """
