rule evo2_run_vep:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/genome.fa.gz",
    output:
        "results/dataset/{dataset}/features/evo2_{version}.parquet",
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
