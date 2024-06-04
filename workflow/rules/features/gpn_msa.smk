rule gpn_msa_run_vep_llr:
    output:
        "results/features/{dataset}/GPN-MSA_LLR.parquet",
    threads:
        workflow.cores
    shell:
        """
        torchrun --nproc_per_node $(echo $CUDA_VISIBLE_DEVICES | awk -F',' '{{print NF}}') \
        -m gpn.msa.inference vep {wildcards.dataset} {config[gpn_msa][msa_path]} \
        {config[gpn_msa][window_size]} {config[gpn_msa][model_path]} {output} \
        --per_device_batch_size 2048 --dataloader_num_workers {threads}
        """


rule gpn_msa_abs_llr:
    input:
        "results/features/{dataset}/GPN-MSA_LLR.parquet",
    output:
        "results/features/{dataset}/GPN-MSA_absLLR.parquet",
    run:
        df = pd.read_parquet(input[0])
        df = df.abs()
        df.to_parquet(output[0], index=False)


rule gpn_msa_run_vep_inner_products:
    output:
        "results/features/{dataset}/GPN-MSA_InnerProducts.parquet",
    threads:
        workflow.cores
    shell:
        """
        torchrun --nproc_per_node $(echo $CUDA_VISIBLE_DEVICES | awk -F',' '{{print NF}}') \
        -m gpn.msa.inference vep_embedding {wildcards.dataset} {config[gpn_msa][msa_path]} \
        {config[gpn_msa][window_size]} {config[gpn_msa][model_path]} {output} \
        --per_device_batch_size 2048 --dataloader_num_workers {threads}
        """