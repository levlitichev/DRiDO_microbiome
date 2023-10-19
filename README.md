# DRiDO Microbiome Study

_Author_: Lev Litichevskiy  
_Last updated_: October 2023  

bioRxiv: XXX

This repository contains code and data related to the DRiDO microbiome manuscript. The starting point for this repository is summarized tables of taxonomic and functional classification results, not fastq files.

## Data processing prior to this repository

Quality control:

- Sunbeam, including cutadapt to remove adapters, trimmomatic for quality-trimming, komplexity to remove highly repetitive sequences, and bowtie2 to remove host reads

Taxonomic classification:

- Kraken2 + MGBC
- HUMAnN3 (MetaPhlAn4)

Functional classification:

- HUMAnN3

## Table of contents

- `data/`
    - `metadata/` 
- `scripts/`: contains data processing scripts, a mix of .R and .Rmd files
- `analysis/`: contains R notebooks used for making figures

## Quick start

The most useful data tables are probably the following:

- `kraken_matrix_agg_by_stool_ID_n1303x2997.txt` with `kraken_taxonomy_n1303.txt`
- `pathabundance_tpm_agg_by_stool_ID_n422x2997.txt`

## Details

Metadata stored separately. Note that sequencing metadata, library metadata, and stool metadata is stored separately. Multiple sequencing IDs (`seq.ID`) can correspond to the same library ID (`lib.ID`), and multiple library IDs can correspond to the same stool sample (`stool.ID`). Each mouse contributed multiple stool samples.


