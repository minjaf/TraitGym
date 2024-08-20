rule gpn_msa_run_vep_llr:
    input:
        "results/dataset/{dataset}/test.parquet",
    output:
        "results/dataset/{dataset}/features/GPN-MSA_LLR.parquet",
    threads:
        workflow.cores
    shell:
        """
        python \
        -m gpn.msa.inference vep {input} {config[gpn_msa][msa_path]} \
        {config[gpn_msa][window_size]} {config[gpn_msa][model_path]} {output} \
        --per_device_batch_size 2048 --dataloader_num_workers {threads} --is_file
        """


rule gpn_msa_abs_llr:
    input:
        "results/dataset/{dataset}/features/GPN-MSA_LLR.parquet",
    output:
        "results/dataset/{dataset}/features/GPN-MSA_absLLR.parquet",
    run:
        df = pd.read_parquet(input[0])
        df = df.abs()
        df.to_parquet(output[0], index=False)


rule gpn_msa_run_vep_inner_products:
    input:
        "results/dataset/{dataset}/test.parquet",
    output:
        "results/dataset/{dataset}/features/GPN-MSA_InnerProducts.parquet",
    threads:
        workflow.cores
    shell:
        """
        python \
        -m gpn.msa.inference vep_embedding {input} {config[gpn_msa][msa_path]} \
        {config[gpn_msa][window_size]} {config[gpn_msa][model_path]} {output} \
        --per_device_batch_size 2048 --dataloader_num_workers {threads} --is_file
        """


rule gpn_msa_run_vep_influence:
    input:
        "results/dataset/{dataset}/test.parquet",
    output:
        "results/dataset/{dataset}/features/GPN-MSA_Influence.parquet",
    threads:
        workflow.cores
    shell:
        """
        python \
        -m gpn.msa.inference vep_influence {input} {config[gpn_msa][msa_path]} \
        {config[gpn_msa][window_size]} {config[gpn_msa][model_path]} {output} \
        --per_device_batch_size 2048 --dataloader_num_workers {threads} --is_file
        """
