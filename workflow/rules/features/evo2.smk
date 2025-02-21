rule evo2_run_vep_llr:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/genome.fa.gz",
    output:
        "results/dataset/{dataset}/features/Evo2_{version}_LLR.parquet",
    threads:
        workflow.cores
    priority: 20
    shell:
        """
        python workflow/scripts/run_vep_llr_evo2.py {input} \
        {config[evo2][wildcards.version][window_size]} \
        {config[evo2][wildcards.version][model_path]} {output} \
        --is_file --dataloader_num_workers 8 \
        --per_device_batch_size \
        {config[evo2][wildcards.version][per_device_batch_size]} \
        """


rule evo2_run_vep_embeddings:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/genome.fa.gz",
    output:
        "results/dataset/{dataset}/features/Evo2_{version}_LLR.parquet",
    threads:
        workflow.cores
    priority: 20
    shell:
        """
        python workflow/scripts/run_vep_embeddings_evo2.py {input} \
        {config[evo2][wildcards.version][window_size]} \
        {config[evo2][wildcards.version][model_path]} {output} \
        --is_file --dataloader_num_workers 8 \
        --per_device_batch_size \
        {config[evo2][wildcards.version][per_device_batch_size]} \
        """
