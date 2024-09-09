#"--per-device-batch-size 32",


#rule run_vep_nucleotide_transformer:
#    input:
#        "results/genome.fa.gz",
#    output:
#        "results/preds/{dataset}/{model}.parquet",
#    wildcard_constraints:
#        dataset="|".join(datasets + ["results/variants_enformer", "results/gnomad/all/defined/128"]),
#        model="|".join(nucleotide_transformer_models)
#    threads:
#        workflow.cores
#    params:
#        lambda wildcards: nucleotide_transformer_params[wildcards.model],
#    priority: 20
#    shell:
#        """
#        python workflow/scripts/run_vep_nucleotide_transformer.py {wildcards.dataset} {input} \
#        {wildcards.model} {output} --dataloader-num-workers 8 {params}
#        """

# seems to be running out of memory, will not spend more time on this since it's not a priority
# torchrun --nproc_per_node=$(echo $CUDA_VISIBLE_DEVICES | awk -F',' '{{print NF}}') workflow/scripts/run_vep_nucleotide_transformer.py {wildcards.dataset} {input} \

rule nucleotide_transformer_run_vep_inner_products:
    input:
        "results/dataset/{dataset}/test.parquet",
        "results/genome.fa.gz",
    output:
        "results/dataset/{dataset}/features/NucleotideTransformer_InnerProducts.parquet",
    threads:
        workflow.cores
    priority: 20
    shell:
        """
        python workflow/scripts/run_vep_embeddings_nucleotide_transformer.py \
        {input} {config[nucleotide_transformer][model_path]} {output} \
        --is_file --dataloader_num_workers 8 --per_device_batch_size 32
        """
