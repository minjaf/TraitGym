rule gpn_run_vep_llr:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/genome.fa.gz",
    output:
        "results/dataset/{dataset}/features/GPN_LLR.parquet",
    threads:
        workflow.cores
    priority: 102
    shell:
        """
        python \
        -m gpn.ss.run_vep {input} {config[gpn][window_size]} \
        {config[gpn][model_path]} {output} --is_file \
        --per_device_batch_size 2048 --dataloader_num_workers {threads}
        """


rule gpn_run_vep_inner_products:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/genome.fa.gz",
    output:
        "results/dataset/{dataset}/features/GPN_InnerProducts.parquet",
    threads:
        workflow.cores
    priority: 102
    shell:
        """
        python \
        -m gpn.ss.run_vep_inner_products {input} {config[gpn][window_size]} \
        {config[gpn][model_path]} {output} --is_file \
        --per_device_batch_size 2048 --dataloader_num_workers {threads}
        """


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


rule gpn_version_run_vep_llr:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/genome.fa.gz",
        lambda wildcards: config["gpn"][wildcards.version]["model_path"],
    output:
        "results/dataset/{dataset}/features/GPN_{version}_LLR.parquet",
    threads:
        workflow.cores
    priority: 102
    shell:
        """
        python \
        -m gpn.ss.run_vep {input[0]} {input[1]} {config[gpn][window_size]} \
        {input[2]} {output} --is_file \
        --per_device_batch_size {config[gpn][per_device_batch_size]} \
        --dataloader_num_workers {threads}
        """


rule gpn_version_run_inner_products:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/genome.fa.gz",
        lambda wildcards: config["gpn"][wildcards.version]["model_path"],
    output:
        "results/dataset/{dataset}/features/GPN_{version}_InnerProducts.parquet",
    threads:
        workflow.cores
    priority: 102
    shell:
        """
        python \
        -m gpn.ss.run_vep_inner_products {input[0]} {input[1]} {config[gpn][window_size]} \
        {input[2]} {output} --is_file \
        --per_device_batch_size {config[gpn][per_device_batch_size]} \
        --dataloader_num_workers {threads}
        """