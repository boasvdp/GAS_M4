# **Streptococcus pyogenes** M4

## Introduction

This repo describes the workflow used to analyse S. pyogenes emm4 genomes from the Netherlands and other countries. The pipeline stored here is meant for this specific project.

Publication pending

## Methods

To use this pipeline, install `snakemake` and `conda` (or better yet, `mamba`).

Download our data from ENA and store the read files in the `reads` directory.

If you want to include public surveillance data from the USA and Norway as well, run the bash script to download this data (updated Dec2022):

```bash
bash workflow/scripts/download_public_data.sh
```

Then, run the snakemake pipeline when snakemake and conda/mamba are installed:

```bash
# Run assembly pipeline
snakemake -s workflow/assembly.smk -j 16 --use-conda

# Perform QC on assembly results
python workflow/scripts/qc.py --species "Streptococcus pyogenes" --output "backup/qc_report_spyo.tsv"

# Run typing pipeline
snakemake -s workflow/typing.smk -j 16 --use-conda
```
