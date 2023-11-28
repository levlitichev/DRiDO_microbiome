# DRiDO Microbiome Study

_Author_: Lev Litichevskiy  
_Last updated_: November 27, 2023  

This repository contains code and data that can be used to reproduce figures from the Dietary Restriction in Diversity Outbred mice (DRiDO) microbiome manuscript. The starting point for this repository is summarized tables of taxonomic and functional classification results (not fastq files).

## Pre-processing of DO metagenomic sequencing data

We used [Sunbeam](https://github.com/sunbeam-labs/sunbeam) for quality control of sequencing reads. More specifically, we used cutadapt to remove adapters, trimmomatic for quality-trimming, [komplexity](https://github.com/eclarke/komplexity) to remove highly repetitive sequences, and bowtie2 to remove host (mm10) reads. We performed taxonomic classification with Kraken2 against the  index available from the [Mouse Gastrointestinal Bacterial Catalogue (MGBC)](https://github.com/BenBeresfordJones/MGBC), and we performed functional classification using HUMAnN3.

## Repository organization

* `analysis` contains .Rmd notebooks used for generating figures
* `scripts` contains a mix of .R and .Rmd files used for data processing and running linear models
* `results` contains the output of `scripts`
* `data` contains the inputs to `scripts` and `analysis`, including metadata

Note that there are multiple layers of metadata: sequencing metadata, library metadata, and stool metadata are stored separately. Multiple sequencing IDs (`seq.ID`) can correspond to the same library ID (`lib.ID`), and multiple library IDs can correspond to the same stool sample (`stool.ID`). Mice contributed one or more stool samples.

## Quick start

### Taxonomy and pathway tables

- Taxonomy: `data/kraken_matrix_agg_by_stool_ID_n1303x2997.txt` with `data/kraken_taxonomy_n1303.txt`
  - Absolute counts
- Pathways: `data/pathabundance_tpm_agg_by_stool_ID_n422x2997.txt`
  - TPM (transcripts-per-million) abundances
  - That is, normalized for gene length and sequencing depth

### Tutorial

This [tutorial](analysis/tutorial.md) demonstrates how to import taxonomic data and perform several basic analyses.
