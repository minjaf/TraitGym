# TraitGym

## Quick start
- Load a dataset
    ```python
    from datasets import load_dataset
    
    dataset = load_dataset("songlab/TraitGym", "mendelian_traits", split="test")
    ```
- Example notebook to run variant effect prediction with a gLM, runs in 5 min on Google Colab: `TraitGym.ipynb` [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/songlab-cal/TraitGym/blob/main/TraitGym.ipynb)

## Resources (https://huggingface.co/datasets/songlab/TraitGym)

## Code (https://github.com/songlab-cal/TraitGym)
- Tries to follow [recommended Snakemake structure](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html)
