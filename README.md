# 🧬 TraitGym
[Benchmarking DNA Sequence Models for Causal Regulatory Variant Prediction in Human Genetics](https://www.biorxiv.org/content/10.1101/2025.02.11.637758v1)

🏆 Leaderboard: https://huggingface.co/spaces/songlab/TraitGym-leaderboard

## ⚡️ Quick start
- Load a dataset
    ```python
    from datasets import load_dataset
    
    dataset = load_dataset("songlab/TraitGym", "mendelian_traits", split="test")
    ```
- Example notebook to run variant effect prediction with a gLM, runs in 5 min on Google Colab: `TraitGym.ipynb` [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/songlab-cal/TraitGym/blob/main/TraitGym.ipynb)

## 🤗 Resources (https://huggingface.co/datasets/songlab/TraitGym)
- Datasets: `{dataset}/test.parquet`
- Subsets: `{dataset}/subset/{subset}.parquet`
- Features: `{dataset}/features/{features}.parquet`
- Predictions: `{dataset}/preds/{subset}/{model}.parquet`
- Metrics: `{dataset}/{metric}/{subset}/{model}.csv`

`dataset` examples (`load_dataset` config name):
- `mendelian_traits_matched_9` (`mendelian_traits`)
- `complex_traits_matched_9` (`complex_traits`)
- `mendelian_traits_all` (`mendelian_traits_full`)
- `complex_traits_all` (`complex_traits_full`)

`subset` examples:
- `all` (default)
- `3_prime_UTR_variant`
- `disease`
- `BMI`

`features` examples:
- `GPN-MSA_LLR`
- `GPN-MSA_InnerProducts`
- `Borzoi_L2`

`model` examples:
-  `GPN-MSA_LLR.minus.score`
-  `GPN-MSA.LogisticRegression.chrom`
-  `CADD+GPN-MSA+Borzoi.LogisticRegression.chrom`

`metric` examples:
- `AUPRC_by_chrom_weighted_average` (main metric)
- `AUPRC`

## 💻 Code (https://github.com/songlab-cal/TraitGym)
- Tries to follow [recommended Snakemake structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)
- GPN-Promoter code is in [the main GPN repo](https://github.com/songlab-cal/gpn)

### Installation
First, clone the repo and `cd` into it.  
Second, install the dependencies:
```bash
conda env create -f workflow/envs/general.yaml
conda activate TraitGym
```
Optionally, download precomputed datasets and predictions (6.7G):
```bash
mkdir -p results/dataset
huggingface-cli download songlab/TraitGym --repo-type dataset --local-dir results/dataset/
```

### Running
To compute a specific result, specify its path:
```bash
snakemake --cores all <path>
```
Example paths (these are already computed):
```bash
# zero-shot LLR
results/dataset/complex_traits_matched_9/AUPRC_by_chrom_weighted_average/all/GPN-MSA_absLLR.plus.score.csv
# logistic regression/linear probing
results/dataset/complex_traits_matched_9/AUPRC_by_chrom_weighted_average/all/GPN-MSA.LogisticRegression.chrom.csv
```
We recommend the following:
```bash
# Snakemake sometimes gets confused about which files it needs to rerun and this forces
# not to rerun any existing file
snakemake --cores all <path> --touch
# to output an execution plan
snakemake --cores all <path> --dry-run
```
To evaluate your own set of model features, place a dataframe of shape `n_variants,n_features` in `results/dataset/{dataset}/features/{features}.parquet`.  
For zero-shot evaluation of column `{feature}` and sign `{sign}` (`plus` or `minus`), you would invoke:
```bash
snakemake --cores all results/dataset/{dataset}/{metric}/all/{features}.{sign}.{feature}.csv
```
To train and evaluate a logistic regression model, you would invoke:
```bash
snakemake --cores all results/dataset/{dataset}/{metric}/all/{feature_set}.LogisticRegression.chrom.csv
```
where `{feature_set}` should first be defined in `feature_sets` in `config/config.yaml` (this allows combining features defined in different files).

## Citation
[Link to paper](https://www.biorxiv.org/content/10.1101/2025.02.11.637758v1)
```bibtex
@article{traitgym,
	author = {Benegas, Gonzalo and Eraslan, Gokcen and Song, Yun S.},
	title = {Benchmarking DNA Sequence Models for Causal Regulatory Variant Prediction in Human Genetics},
	elocation-id = {2025.02.11.637758},
	year = {2025},
	doi = {10.1101/2025.02.11.637758},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2025/02/12/2025.02.11.637758},
	eprint = {https://www.biorxiv.org/content/early/2025/02/12/2025.02.11.637758.full.pdf},
	journal = {bioRxiv}
}
```
