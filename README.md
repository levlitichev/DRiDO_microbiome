DRiDO Microbiome Study
================

*Author*: Lev Litichevskiy  
*Last updated*: December 18, 2024

This repository contains code and data that can be used to reproduce
figures from the Dietary Restriction in Diversity Outbred mice (DRiDO)
microbiome manuscript. The starting point for this repository is
summarized tables of taxonomic and functional classification results
(not fastq files).

## 1. Citation

L Litichevskiy, M Considine, J Gill, V Shandar, … A Di Francesco, GA Churchill, M Li, CA Thaiss. Interactions between the gut microbiome, dietary restriction, and aging in genetically diverse mice. <https://www.biorxiv.org/content/10.1101/2023.11.28.568137>

## 2. Quick start

### Data

- Taxonomy: `data/kraken_matrix_agg_by_stool_ID_n1303x2997.txt` with
  `data/kraken_taxonomy_n1303.txt`
  - Absolute counts
- Pathways: `data/pathabundance_tpm_agg_by_stool_ID_n422x2997.txt`
  - TPM (transcripts-per-million) abundances
  - That is, normalized for gene length and sequencing depth
- fastq files available on SRA: [PRJNA1054518](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1054518)

### Tutorial

This [tutorial](analysis/tutorial.md) demonstrates how to import
taxonomic data and perform several basic analyses.

## 3. Repository organization

- `analysis` contains .Rmd notebooks used for generating figures
- `plots` contains figures, i.e. the output of `analysis`
- `scripts` contains a mix of .R and .Rmd files used for data processing
  and running linear models
- `results` contains the output of `scripts`
- `data` contains the inputs to `scripts` and `analysis`, including
  metadata

See [here](TOC.md) for more details about the overall workflow.

See [here](script_used_for_each_figure_panel.md) for which script was used to produce every figure panel in the manuscript.

Note that there are multiple layers of metadata: sequencing metadata,
library metadata, and stool metadata are stored separately. Multiple
sequencing IDs (`seq.ID`) can correspond to the same library ID
(`lib.ID`), and multiple library IDs can correspond to the same stool
sample (`stool.ID`). Mice contributed one or more stool samples.

Note also that library and stool metadata is embedded in the SRA metadata, but sequencing metadata is not. This is because individual SRA runs correspond to (unique) library IDs. In the unlikely event that a user will need to know which sequencing IDs correspond to each SRA run (i.e. SRR accession), we have made available a [mapping](https://github.com/levlitichev/DRiDO_microbiome/blob/main/data/metadata/library_metadata_after_QC_n3213_230620.txt) between SRR accessions and sequencing IDs. Only a small number of libraries / SRA runs (n=259) correspond to multiple sequencing IDs, i.e. these libraries were sequenced multiple times.

## 4. System requirements

This code was run on macOS Big Sur using R v4.2.2. All R packages are
available from CRAN or Bioconductor — except for ASReml, which requires
a license. ASReml was used for estimating heritability and running
linear mixed models. Identical results can be produced using the lme4qtl
package (see [run_lme4qtl.R](scripts/run_lme4qtl.R) for an example).

All analyses except for mediation and QTL mapping were run on a laptop.
Mediation analysis was performed on a cluster using
[Snakemake](https://snakemake.github.io/)
([Snakefile_mediation](scripts/Snakefile_mediation),
[run_mediation_one_diet_one_pheno.R](scripts/run_mediation_one_diet_one_pheno.R)),
and QTL mapping was performed on a cluster using R/qtl2
([run_genetic_mapping_rqtl2.R](scripts/run_genetic_mapping_rqtl2.R)).

## 5. QTL mapping

Karl Broman’s [R/qtl2](https://kbroman.org/qtl2/) was written
specifically to handle multi-parent QTL mapping crosses such as DO mice.

QTL mapping was performed as described in [Zhang et al., *Genetics*,
2022](https://academic.oup.com/genetics/article/220/1/iyab157/6375446#325918956)
(“Genetic linkage analysis”).

Input files not included in this repo:

1.  `cc_variants.sqlite`: Imputed variants from DO founders. Download
    [here](https://figshare.com/articles/dataset/SQLite_database_of_variants_in_Collaborative_Cross_founder_mouse_strains/5280229).

2.  `prob.8state.allele.qtl2_200131.Rdata`: 8 state allele probabilities
    for DRiDO mice. Download
    [here](https://figshare.com/articles/dataset/Supplementary_files_associated_with_the_DRiDO_microbiome_manuscript/25043753).
