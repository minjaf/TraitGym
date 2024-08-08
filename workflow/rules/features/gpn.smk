rule gpn_run_vep_embed_dists:
    input:
        "results/genome.fa.gz",
    output:
        "results/features/{dataset}/GPN_EmbedDists.parquet",
    threads:
        workflow.cores
    shell:
        """
        python \
        -m gpn.ss.run_vep_embed_dists {wildcards.dataset} {input} \
        {config[gpn][window_size]} {config[gpn][model_path]} {output} \
        --per-device-batch-size 2048 --dataloader-num-workers {threads}
        """
