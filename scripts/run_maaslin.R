# Run linear mixed models (LMMs) with MaAsLin.
# Already prepared matrices of features that I want to test, including community features.
# No prevalence filter, no transformation, no normalization.
# The only thing I need to do before running the LMMs is scale all features
# so the coefficients for the community features are on the same as those
# for the individual features.

# load libraries
library(tidyverse)
library(Maaslin2)

# change working directory
setwd("~/DRiDO_microbiome_github/scripts")

# --- import metadata ---

# DO
DO.stool.meta.df <- read.table(
  "../data/metadata/stool_metadata_after_QC_no_controls_n2997_240418.txt", 
  sep="\t", header=T)
DO.mouse.meta.df <- read.csv(
  "../data/metadata/AnimalData_Processed_20230712.csv")
DO.meta.df <- DO.stool.meta.df %>% 
  merge(DO.mouse.meta.df, by.x="mouse.ID", by.y="MouseID")

# subset to AL
DO.AL.meta.df <- DO.meta.df %>% 
  
  dplyr::filter(Diet == "AL") %>% 
  
  # sort by stool.ID, then set as rownames
  arrange(stool.ID) %>%
  column_to_rownames("stool.ID") %>%

  # rename so the names show up nicely in the model output
  mutate(Age = age.approx.months) %>% 
  mutate(Cage = paste0("c", HID)) %>% 
  mutate(Batch = ext.batch)

# B6
B6.meta.df <- read.table(
  "../data/metadata/B6_sample_metadata.txt",
  sep="\t", header=T) %>% 
  
  # rename so the names show up nicely in the model output
  dplyr::rename(Age = age.months, Cage = cage) %>% 
  
  # sort by id, then set as rownames
  arrange(id) %>% 
  column_to_rownames("id")

# Human
human.meta.df <- read.table(
  "../data/metadata/CMD_human_sample_metadata_n4101.txt",
  sep="\t", header=T) %>% 
  
  # sort by sample_id, then set as rownames
  arrange(sample_id) %>% 
  column_to_rownames("sample_id")

# --- DO AL, genera ---

# import microbiome data
DO.AL.genus.df.feats.in.rows <- read.table(
  "../results/DO_AL_genus_log2relab_filt_w_comm_n253x573.txt",
  sep="\t", header=T, row.names=1)
dim(DO.AL.genus.df.feats.in.rows)

# transpose and scale features
DO.AL.genus.df <- DO.AL.genus.df.feats.in.rows %>% 
  
  # transpose because MaAsLin2 wants samples in the rows
  t() %>%
  
  # scale
  scale() %>% 
  
  # convert back to dataframe
  data.frame(check.names=F) %>% 
  
  # sort samples alphabetically, then set back to rownames
  rownames_to_column("stool.ID") %>% 
  arrange(stool.ID) %>% 
  column_to_rownames("stool.ID")

# run MaAsLin2
maaslin.res.DO.AL.genus <- Maaslin2(
  input_data = DO.AL.genus.df,
  input_metadata = DO.AL.meta.df,
  output = "../results/maaslin2_DO_AL_genus",
  fixed_effects = c("Age"),
  random_effects = c("mouse.ID", "Cohort", "Cage", "Batch"),
  cores = 4,
  min_abundance = 0, min_prevalence = 0,
  normalization = "NONE", transform="NONE",
  plot_heatmap = FALSE,
  plot_scatter = FALSE)

# --- DO AL, pathways ---

# import microbiome data
DO.AL.pathway.df.feats.in.rows <- read.table(
  "../results/DO_AL_pathway_log2tpm_filt_w_comm_n263x573.txt",
  sep="\t", header=T, quote="", row.names=1)
dim(DO.AL.pathway.df.feats.in.rows)

# transpose and scale features
DO.AL.pathway.df <- DO.AL.pathway.df.feats.in.rows %>% 
  
  # transpose because MaAsLin2 wants samples in the rows
  t() %>%
  
  # scale
  scale() %>% 
  
  # convert back to dataframe
  data.frame(check.names=F) %>% 
  
  # sort samples alphabetically, then set back to rownames
  rownames_to_column("stool.ID") %>% 
  arrange(stool.ID) %>% 
  column_to_rownames("stool.ID")

# run MaAsLin2
maaslin.res.DO.AL.pathway <- Maaslin2(
  input_data = DO.AL.pathway.df,
  input_metadata = DO.AL.meta.df,
  output = "../results/maaslin2_DO_AL_pathway",
  fixed_effects = c("Age"),
  random_effects = c("mouse.ID", "Cohort", "Cage", "Batch"),
  cores = 4,
  min_abundance = 0, min_prevalence = 0,
  normalization = "NONE", transform="NONE",
  plot_heatmap = FALSE,
  plot_scatter = FALSE)

# --- B6, genera ---

# import microbiome data
B6.genus.df.feats.in.rows <- read.table(
  "../results/B6_genus_log2relab_filt_w_comm_n267x141.txt",
  sep="\t", check.names=F, header=T, row.names=1)
dim(B6.genus.df.feats.in.rows)

# transpose and scale features
B6.genus.df <- B6.genus.df.feats.in.rows %>% 
  
  # transpose because MaAsLin2 wants samples in the rows
  t() %>%
  
  # scale
  scale() %>% 
  
  # convert back to dataframe
  data.frame(check.names=F) %>% 
  
  # sort samples alphabetically, then set back to rownames
  rownames_to_column("id") %>% 
  arrange(id) %>% 
  column_to_rownames("id")


# run MaAsLin2
maaslin.res.B6.genus <- Maaslin2(
  input_data = B6.genus.df,
  input_metadata = B6.meta.df,
  output = "../results/maaslin2_B6_genus",
  fixed_effects = c("Age"),
  random_effects = c("Cage"),
  cores = 4,
  min_abundance = 0, min_prevalence = 0,
  normalization = "NONE", transform="NONE",
  plot_heatmap = FALSE,
  plot_scatter = FALSE)

# --- B6, pathways ---

# import microbiome data
B6.pathway.df.feats.in.rows <- read.table(
  "../results/B6_pathway_log2tpm_filt_w_comm_n234x141.txt",
  sep="\t", check.names=F, header=T, quote="", row.names=1)
dim(B6.pathway.df.feats.in.rows)

# transpose and scale features
B6.pathway.df <- B6.pathway.df.feats.in.rows %>% 
  
  # transpose because MaAsLin2 wants samples in the rows
  t() %>%
  
  # scale
  scale() %>% 
  
  # convert back to dataframe
  data.frame(check.names=F) %>% 
  
  # sort samples alphabetically, then set back to rownames
  rownames_to_column("id") %>% 
  arrange(id) %>% 
  column_to_rownames("id")

# run MaAsLin2
maaslin.res.B6.pathway <- Maaslin2(
  input_data = B6.pathway.df,
  input_metadata = B6.meta.df,
  output = "../results/maaslin2_B6_pathway",
  fixed_effects = c("Age"),
  random_effects = c("Cage"),
  cores = 4,
  min_abundance = 0, min_prevalence = 0,
  normalization = "NONE", transform="NONE",
  plot_heatmap = FALSE,
  plot_scatter = FALSE)

# --- Human, genera ---

# import microbiome data
human.genus.df.feats.in.rows <- read.table(
  "../results/CMD_human_genus_log2relab_filt_w_comm_n95x4101.txt",
  sep="\t", check.names=F, header=T, row.names=1)

# transpose and scale features
human.genus.df <- human.genus.df.feats.in.rows %>% 
  
  # transpose because MaAsLin2 wants samples in the rows
  t() %>%
  
  # scale
  scale() %>% 
  
  # convert back to dataframe
  data.frame(check.names=F) %>% 
  
  # sort samples alphabetically, then set back to rownames
  rownames_to_column("sample_id") %>% 
  arrange(sample_id) %>% 
  column_to_rownames("sample_id")

# run MaAsLin2
maaslin.res.human.genus <- Maaslin2(
  input_data = human.genus.df,
  input_metadata = human.meta.df,
  output = "../results/maaslin2_human_genus",
  fixed_effects = c("age"),
  random_effects = c("study_name"),
  cores = 4,
  min_abundance = 0, min_prevalence = 0,
  normalization = "NONE", transform="NONE",
  plot_heatmap = FALSE,
  plot_scatter = FALSE)


# --- Human, pathways ---

# import microbiome data
human.pathway.df.feats.in.rows <- read.table(
  "../results/CMD_human_pathway_log2tpm_filt_w_comm_n359x4101.txt",
  sep="\t", check.names=F, quote="", header=T, row.names=1)

# transpose and scale features
human.pathway.df <- human.pathway.df.feats.in.rows %>% 
  
  # transpose because MaAsLin2 wants samples in the rows
  t() %>%
  
  # scale
  scale() %>% 
  
  # convert back to dataframe
  data.frame(check.names=F) %>% 
  
  # sort samples alphabetically, then set back to rownames
  rownames_to_column("sample_id") %>% 
  arrange(sample_id) %>% 
  column_to_rownames("sample_id")

# run MaAsLin2
maaslin.res.human.pathway <- Maaslin2(
  input_data = human.pathway.df,
  input_metadata = human.meta.df,
  output = "../results/maaslin2_human_pathway",
  fixed_effects = c("age"),
  random_effects = c("study_name"),
  cores = 4,
  min_abundance = 0, min_prevalence = 0,
  normalization = "NONE", transform="NONE",
  plot_heatmap = FALSE,
  plot_scatter = FALSE)
