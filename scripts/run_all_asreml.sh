# This script will run all linear mixed models (LMMs) -- except those run with
# MaAsLin (see maaslin.Rmd) -- that evaluate how different experimental
# variables affect microbiome features. This script does not
# run the microbiome-phenotype associations.

cd ~/DRiDO_microbiome_github/scripts

# --- TAXONOMY ---

# Model 1: Kraken, model with time, genera, fixed and random effects
# y_mb ~ age (f) + DR (f) + time (r) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
in_path="../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt"
out_dir="../results/asreml_kraken_genus_w_time_ranef/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Quarter.Date + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Quarter.Date + Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

# Model 2: Model with time as a fixed effect
# y_mb ~ age (f) + DR (f) + time (f) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
in_path="../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt"
out_dir="../results/asreml_kraken_genus_w_time_fixef/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ age.wks.scaled + Diet.5mo.as.AL + time.scaled" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

# Model 3: Model without time
# y_mb ~ age (f) + DR (f) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
in_path="../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt"
out_dir="../results/asreml_kraken_genus_wo_time/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

# Model 4: Run separately for each of three cross-sectional data slices
# y_mb ~ age (f) + DR (f) + [genetics] (r) + batch (r) + cohort (r) + cage (r)

## Slice 1
in_path="../results/cross_sectional_slices/kraken_genus_clr_filt_w_comm_slice1_n107x573.txt"
out_dir="../results/cross_sectional_slices/asreml_kraken_genus_slice1/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

## Slice 2
in_path="../results/cross_sectional_slices/kraken_genus_clr_filt_w_comm_slice2_n107x424.txt"
out_dir="../results/cross_sectional_slices/asreml_kraken_genus_slice2/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

## Slice 3
in_path="../results/cross_sectional_slices/kraken_genus_clr_filt_w_comm_slice3_n107x381.txt"
out_dir="../results/cross_sectional_slices/asreml_kraken_genus_slice3/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

# Model 9: Model run separately per age
# At 5 months:          y_mb ~ [genetics] (r) + batch (r) + cohort (r) + cage (r) 
# After randomization:  y_mb ~ DR (f) + [genetics] (r) + batch (r) + cohort (r) + cage (r) 

## 5 mo
in_path="../results/split_by_age/kraken_genus_clr_filt_w_comm_5mo_n107x569.txt"
out_dir="../results/split_by_age/asreml_kraken_genus_5mo/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ 1" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

## 10 mo
in_path="../results/split_by_age/kraken_genus_clr_filt_w_comm_10mo_n107x513.txt"
out_dir="../results/split_by_age/asreml_kraken_genus_10mo/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

## 16 mo
in_path="../results/split_by_age/kraken_genus_clr_filt_w_comm_16mo_n107x646.txt"
out_dir="../results/split_by_age/asreml_kraken_genus_16mo/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

## 22 mo
in_path="../results/split_by_age/kraken_genus_clr_filt_w_comm_22mo_n107x522.txt"
out_dir="../results/split_by_age/asreml_kraken_genus_22mo/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

## 28 mo
in_path="../results/split_by_age/kraken_genus_clr_filt_w_comm_28mo_n107x368.txt"
out_dir="../results/split_by_age/asreml_kraken_genus_28mo/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + Cohort + Batch + Cage" \
  "~ Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir # also deletes the uncollated .txt files
rm "$out_dir"*.Rds

# Model 10: All random effects
in_path="../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt"
out_dir="../results/asreml_kraken_genus_all_ranef_w_time/"
mkdir $out_dir
Rscript run_asreml_all_ranef_w_time.R $in_path $out_dir
Rscript collate_asreml_all_ranef.R $out_dir
find "$out_dir" -type f -name "*.txt" ! -name "ranef.txt" -delete

# Model 1, with MetaPhlAn instead of Kraken taxonomic results
# y_mb ~ age (f) + DR (f) + time (r) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
in_path="../results/metaphlan_genus_log2relab_filt_w_comm_n259x2997.txt"
out_dir="../results/asreml_metaphlan_genus_w_time_ranef/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Quarter.Date + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Quarter.Date + Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

# Model 1, with Kraken species-level data
# y_mb ~ age (f) + DR (f) + time (r) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
in_path="../results/kraken_species_clr_filt_w_comm_n207x2997.txt"
out_dir="../results/asreml_kraken_species_w_time_ranef/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Quarter.Date + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Quarter.Date + Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

# Model 1, with downsampled data
# y_mb ~ age (f) + DR (f) + time (r) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
in_path="../results/kraken_genus_clr_filt_w_comm_downsampled_n107x519.txt"
out_dir="../results/asreml_kraken_genus_downsampled/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Quarter.Date + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Quarter.Date + Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

# --- PATHWAYS ---

# Model 1: HUMAnN pathways, model with time, fixed and random effects
# y_mb ~ age (f) + DR (f) + time (r) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r)
in_path="../results/pathway_log2tpm_filt_w_comm_n273x2997.txt"
out_dir="../results/asreml_humann_pathways_w_time_ranef/"
mkdir $out_dir
Rscript run_asreml.R $in_path $out_dir \
  "~ age.wks.scaled + Diet.5mo.as.AL" \
  "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Quarter.Date + Cohort + Cage + Batch" \
  "~ ide(Mouse, kinship.mat.x2) + Quarter.Date + Cohort + Cage + Batch"
Rscript collate_asreml.R $out_dir
rm "$out_dir"*.Rds

# Model 10: All random effects
in_path="../results/pathway_log2tpm_filt_w_comm_n273x2997.txt"
out_dir="../results/asreml_humann_pathways_all_ranef_w_time/"
mkdir $out_dir
Rscript run_asreml_all_ranef_w_time.R $in_path $out_dir
Rscript collate_asreml_all_ranef.R $out_dir
find "$out_dir" -type f -name "*.txt" ! -name "ranef.txt" -delete
