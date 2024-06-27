# This script will run all linear mixed models (LMMs) that evaluate how different
# experimental variables affect microbiome features. This script does not
# run the microbiome-phenotype associations.

setwd("~/DRiDO_microbiome_github/scripts/")

# --- TAXONOMY ---

# Model 1: Kraken, model with time, genera, fixed and random effects
# y_mb ~ age (f) + DR (f) + time (f) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
Rscript run_asreml.R \
  ../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt \
  ../results/asreml_kraken_genus_w_time/ \
  "~ age.wks.scaled + Diet.5mo.as.AL + time.scaled" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/asreml_kraken_genus_w_time/

# Model 2: Model without time
# y_mb ~ age (f) + DR (f) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
Rscript run_asreml.R \
  ../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt \
  ../results/asreml_kraken_genus_wo_time/ \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/asreml_kraken_genus_wo_time/

# Model 3: Model run separately for each of three cross-sectional data slices
# y_mb ~ age (f) + DR (f) + [genetics] (r) + batch (r) + cohort (r) + cage (r)

## Slice 1
Rscript run_asreml.R \
  ../results/cross_sectional_slices/kraken_genus_clr_filt_w_comm_slice1_n107x573.txt \
  ../results/cross_sectional_slices/asreml_kraken_genus_slice1/ \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/cross_sectional_slices/asreml_kraken_genus_slice1/

## Slice 2
Rscript run_asreml.R \
  ../results/cross_sectional_slices/kraken_genus_clr_filt_w_comm_slice2_n107x424.txt \
  ../results/cross_sectional_slices/asreml_kraken_genus_slice2/ \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/cross_sectional_slices/asreml_kraken_genus_slice2/

## Slice 3
Rscript run_asreml.R \
  ../results/cross_sectional_slices/kraken_genus_clr_filt_w_comm_slice3_n107x381.txt \
  ../results/cross_sectional_slices/asreml_kraken_genus_slice3/ \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/cross_sectional_slices/asreml_kraken_genus_slice3/

# Model 8: Model run separately per age
# At 5 months:          y_mb ~ [genetics] (r) + batch (r) + cohort (r) + cage (r) 
# After randomization:  y_mb ~ DR (f) + [genetics] (r) + batch (r) + cohort (r) + cage (r) 

## 5 mo
Rscript run_asreml.R \
  ../results/split_by_age/kraken_genus_clr_filt_w_comm_5mo_n107x569.txt \
  ../results/split_by_age/asreml_kraken_genus_5mo/ \
  "~ 1" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/split_by_age/asreml_kraken_genus_5mo/

## 10 mo
Rscript run_asreml.R \
  ../results/split_by_age/kraken_genus_clr_filt_w_comm_10mo_n107x513.txt \
  ../results/split_by_age/asreml_kraken_genus_10mo/ \
  "~ Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/split_by_age/asreml_kraken_genus_10mo/

## 16 mo
Rscript run_asreml.R \
  ../results/split_by_age/kraken_genus_clr_filt_w_comm_16mo_n107x646.txt \
  ../results/split_by_age/asreml_kraken_genus_16mo/ \
  "~ Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/split_by_age/asreml_kraken_genus_16mo/

## 22 mo
Rscript run_asreml.R \
  ../results/split_by_age/kraken_genus_clr_filt_w_comm_22mo_n107x522.txt \
  ../results/split_by_age/asreml_kraken_genus_22mo/ \
  "~ Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/split_by_age/asreml_kraken_genus_22mo/

## 28 mo
Rscript run_asreml.R \
  ../results/split_by_age/kraken_genus_clr_filt_w_comm_28mo_n107x368.txt \
  ../results/split_by_age/asreml_kraken_genus_28mo/ \
  "~ Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/split_by_age/asreml_kraken_genus_28mo/

# Model 9: All random effects
Rscript run_asreml_all_ranef_w_time.R \
  ../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt \
  ../results/asreml_kraken_genus_all_ranef_w_time/
# output is .txt files, not .RDS files, so we don't use collate_asreml.R

# Model 1, with MetaPhlAn instead of Kraken taxonomic results
# y_mb ~ age (f) + DR (f) + time (f) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
Rscript run_asreml.R \
  ../results/metaphlan_genus_log2relab_filt_w_comm_n259x2997.txt \
  ../results/asreml_metaphlan_genus_w_time/ \
  "~ age.wks.scaled + Diet.5mo.as.AL + time.scaled" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/asreml_metaphlan_genus_w_time/

# Model 1, with Kraken species-level data
# y_mb ~ age (f) + DR (f) + time (f) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
Rscript run_asreml.R \
  ../results/kraken_species_clr_filt_w_comm_n207x2997.txt \
  ../results/asreml_kraken_species_w_time/ \
  "~ age.wks.scaled + Diet.5mo.as.AL + time.scaled" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/asreml_kraken_species_w_time/

# Model 1, using lme4qtl instead of ASReml
Rscript run_lme4qtl.R

# Model 1, with downsampled data
# y_mb ~ age (f) + DR (f) + time (f) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
Rscript run_asreml.R \
  ../results/kraken_genus_clr_filt_w_comm_downsampled_n107x519.txt \
  ../results/asreml_kraken_genus_downsampled/ \
  "~ age.wks.scaled + Diet.5mo.as.AL + time.scaled" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
Rscript collate_asreml.R ../results/asreml_kraken_genus_downsampled/

# --- PATHWAYS ---

# Model 1: HUMAnN pathways, model with time, genera, fixed and random effects
run_asreml_w_time.R ../results/DO_pathway_log2tpm_filt_w_comm_n273x2997.txt ../results/asreml_humann_pathways_w_time/

# Model 2: Model without time
run_asreml_wo_time.R ../results/DO_pathway_log2tpm_filt_w_comm_n273x2997.txt ../results/asreml_humann_pathways_wo_time/

# Model 3: Model run separately for samples collected close in chronological time
TODO: run_asreml_at_different_times.R ../results/DO_pathway_log2tpm_filt_w_comm_n273x2997.txt ../results/asreml_humann_pathways_per_time_slice/ # TODO

# Model 8: Model run separately per age group
# probably can be consolidated to use run_asreml_w_time.R and a different data file
run_asreml_per_age.R ../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt ../results/asreml_humann_pathways_per_age_w_time

# Model 9: All random effects
run_asreml_all_ranef_w_time.R ../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt ../results/asreml_kraken_genus_all_ranef_w_time/

# Model 1, with MetaPhlAn instead of Kraken taxonomic results
run_asreml_w_time.R ../results/DO_genus_log2relab_filt_w_comm_n259x2997.txt ../results/asreml_metaphlan_genus_w_time/

# Model 1, with Kraken species-level data
run_asreml_w_time.R results/kraken_species_clr_filt_w_comm_n207x2997_240515.txt ../results/asreml_kraken_species/

# Model 1, using lme4qtl instead of ASReml
run_lme4qtl.R ../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt ../results/lme4qtl_kraken_genus

# Model 1, using downsampled data
# probably can be consolidated to use run_asreml_w_time.R and a different data file
run_asreml_downsampled.R ../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt ../results/asreml_kraken_genus_w_time/