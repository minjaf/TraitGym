# TraitGym
Benchmarking DNA Sequence Models for Causal Regulatory Variant Prediction in Human Genetics

## Quick start
- Load a dataset
    ```python
    from datasets import load_dataset
    
    dataset = load_dataset("songlab/TraitGym", "mendelian_traits", split="test")
    ```
- Example notebook to run variant effect prediction with a gLM, runs in 5 min on Google Colab: `TraitGym.ipynb` [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/songlab-cal/TraitGym/blob/main/TraitGym.ipynb)

## Resources (https://huggingface.co/datasets/songlab/TraitGym)
- Datasets: `{dataset}/test.parquet`
- Subsets: `{dataset}/subset/{subset}.parquet`
- Features: `{dataset}/features/{features}.parquet`
- Predictions: `{dataset}/preds/{subset}/{model}.parquet`
- Metrics: `{dataset}/{metric}/{subset}/{model}.csv`

`{dataset}` examples (`load_dataset` config name):
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

## Code (https://github.com/songlab-cal/TraitGym)
- Tries to follow [recommended Snakemake structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)
