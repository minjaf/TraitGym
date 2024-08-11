rule download_phyloP_100v:
    output:
        "results/conservation/phyloP-100v.bw",
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw -O {output}"


rule download_phastCons_100v:
    output:
        "results/conservation/phastCons-100v.bw",
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.phastCons100way.bw -O {output}"


rule download_phyloP_241m:
    output:
        "results/conservation/phyloP-241m.bw"
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus241way/cactus241way.phyloP.bw -O {output}"


rule conservation_features:
    input:
        "results/conservation/{model}.bw",
    output:
        "results/features/{dataset}/{model,phyloP-100v|phastCons-100v|phyloP-241m}.parquet",
    threads: workflow.cores // 3
    run:
        import pyBigWig

        df = load_dataset(wildcards["dataset"], split="test").to_pandas()
        bw = pyBigWig.open(input[0])
        df["score"] = df.progress_apply(
            lambda v: bw.values(f"chr{v.chrom}", v.pos-1, v.pos)[0],
            axis=1
        )
        df = df[["score"]]
        df.to_parquet(output[0], index=False)
