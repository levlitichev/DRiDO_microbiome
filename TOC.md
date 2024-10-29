### Data preparation (within `scripts`)

| File | Notes|
| ---- | ------ |
| 1. `aggregate_by_stool_ID_kraken.Rmd` | |
| 2. `aggregate_by_stool_ID_humann.Rmd` | |
| 3. `prepare_for_linear_modeling_kraken.Rmd` | |
| 4. `prepare_for_linear_modeling_humann.Rmd` | |
| 5. `extract_human_data_CMD.Rmd` | Needed for `run_maaslin.R`. |
| 6. `create_cross_sectional_slices.Rmd` | Needed to run the linear model at different slices of chronological time (within `run_all_asreml.sh`). |
| 7. `get_phenotype_value_closest_to_microbiome_sample.Rmd` | Needed to perform microbiome-phenotype associations (`run_lme4_mb_pheno_assoc.R`) and mediation analysis (`Snakefile_mediation`). |

### Linear mixed models (within `scripts`)

| File | Notes |
| ---- | ----- |
| `run_all_LMM.sh` | Overview of all LMMs in the paper. |
| →`run_all_asreml.sh` | All LMMs that use ASReml. |
| →→`run_asreml.R` | |
| →→`collate_asreml.R`| |
| →→`run_asreml_all_ranef_w_time.R` | |
| →→`collate_asreml_all_ranef.R` | |
| →`run_maaslin.R` | |
| →`run_lme4qtl.R` | |
| →`collate_lme4qtl.R` | |
| →`run_lme4_mb_pheno_assoc.R` | |
| →`collate_lme4_mb_pheno_assoc.R` | |
| →`lm_cross_sectional_mb_pheno_assoc.Rmd` | |

### Mediation analysis (within `scripts`, must be run on cluster)

| File | Notes |
| ---- | ----- |
| 1. `Snakefile_mediation` | |
| →`run_mediation_one_diet_one_pheno.R` | |
| 2. `aggregate_mediation_results.R` | |

### Analysis and producing plots (within `analysis/`)

N.B. `kraken.Rmd`, `humann.Rmd`, and `prediction.Rmd` account for >50% of the figure panels in the manuscript.

| File | Notes |
| ---- | ----- |
| 1. `qc.Rmd` | |
| 2. `effect_of_chronological_time.Rmd` | |
| 3. `kraken_v_metaphlan.Rmd` | |
| 4. `kraken_unaggregated.Rmd` | |
| 5. `kraken.Rmd` | |
| 6. `metaphlan_poscons.Rmd` | |
| 7. `humann.Rmd` | |
| 8. `permanova.Rmd` | |
| 9. `compare_DO_to_B6_to_humans.Rmd` | |
| 10. `prediction.Rmd` | |
| 11. `cohousing.Rmd` | |
| 12. `perc_heritable_taxa_by_study.Rmd` | |
| 13. `compare_to_schlamp.Rmd` | |
| 14. `asreml_longitudinal_v_cross_sectional.Rmd` | |
| 15. `asreml_v_lme4qtl.Rmd` | |
| 16. `pheno_assoc_and_mediation.Rmd` | | 

### QTL mapping

| File | Notes |
| ---- | ----- |
| 1. `scripts/run_genetic_mapping_rqtl2.R` | Performs QTL mapping, must be done on a cluster. |
| 2. `analysis/QTL.Rmd` | Plots results of QTL mapping. |
