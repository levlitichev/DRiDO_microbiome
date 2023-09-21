# DRiDO Microbiome Study

_Maintainer_: Lev Litichevskiy  
_Last updated_: September 21, 2023  

This repository contains code and data related to the DRiDO microbiome manuscript.

bioRxiv: XXX

## Table of contents

data/
    DO_metaphlan.txt
    DO_pathabundance.txt
    DO_genus_agg_by_stool_ID_nAAxBB.txt
    DO_species_agg_by_stool_ID_nAAxBB.txt
    DO_pathabundance_tpm_agg_by_stool_ID_nAAxBB.txt
    B6/
        TODO
    human/
        TODO       
    metadata/
        TODO
    AllData_20230731.csv
    kinship.all_chroms_downloaded_from_wright22.csv
 
results/
    DO_AL_species_log2relab_filt_w_comm_n292x573.txt
    DO_AL_genus_log2relab_filt_w_comm_n240x573.txt
    DO_AL_pathway_log2tpm_filt_w_comm_n103x573.txt
    DO_species_log2relab_filt_w_comm_nAAxBB.txt
    DO_genus_log2relab_filt_w_comm_nAAxBB.txt
    DO_pathway_log2tpm_filt_w_comm_nAAxBB.txt
    DO_all_features_nAAxBB.txt
    asreml_TODO...
    lme4_mb_pheno_assoc_TODO...
    B6/
        TODO
    human/
        TODO

scripts/
1. qc.Rmd (?) 
2. aggregate_by_stool_ID.Rmd
3. prepare_for_linear_modeling.Rmd
4. run_asreml.R
5. collate_asreml_results.R
6. get_phenotype_value_closest_to_microbiome_sample.Rmd
7. run_lme4_mb_pheno_assoc.R 
8. collate_lme4_mb_pheno_assoc_results.R
9. maaslin.Rmd

analysis/
1. taxonomy_metaphlan.Rmd
2. function_humann.Rmd

