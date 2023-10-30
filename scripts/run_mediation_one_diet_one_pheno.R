# I am using Tingley's "model-based inference" framework. We need to fit two models:
#   M ~ T + X 
#   Y ~ T + M + X
# where Y is a phenotype, M is a microbiome feature, T is diet, and X is covariates.

# I need to run each diet separately (against AL). From the help text within the `mediate` function:
# "The treatment can be either binary (integer or a two-valued factor) or continuous (numeric)."

# Also, I can only provide one random effect. I received this error message when providing two random effects:
#   "mediate does not support more than two levels per model"
# I will provide mouse ID and have to ignore genetics.

# So the models are
#   y_mb ~ diet_X + age + (1|mouse)
#   y_pheno ~ y_mb + diet_X + age + (1|mouse)
#
# This script was run on a cluster using Snakemake.

# Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(mediation))
suppressPackageStartupMessages(library(foreach)) # for parallel for-loop
suppressPackageStartupMessages(library(doParallel)) # for parallel for-loop

# Input
this.diet <- snakemake@params[["diet"]]
this.pheno <- snakemake@params[["pheno"]]
mb.data.path <- snakemake@input[["mb"]]
pheno.data.path <- snakemake@input[["pheno"]]
stool.meta.path <- snakemake@input[["stool_meta"]]
mouse.meta.path <- snakemake@input[["mouse_meta"]]
out.path <- snakemake@output[[1]]
n.threads <- snakemake@threads
log.path <- snakemake@log[[1]]

cat(this.diet, "\n")
cat(this.pheno, "\n")
cat("# threads:", n.threads, "\n")

# Specify mediation and outcome models
med.formula <- as.formula("mb.scaled ~ Diet.5mo.as.AL + age.wks.scaled + (1|mouse.ID)")
out.formula <- as.formula("pheno.scaled ~ mb.scaled + Diet.5mo.as.AL + age.wks.scaled + (1|mouse.ID)")

# Define function that will work on one microbiome feature and one phenotype
run_mediation_one_pheno_one_mb_feat <- function(
    this.pheno, this.mb.feat, this.df, this.diet, med.formula, out.formula) {
  
  # scale microbiome feature
  this.df[["mb.scaled"]] <- as.numeric(scale(this.df[[this.mb.feat]]))
 
  # make sure there are enough unique microbiome values
  if(length(unique(this.df[["mb.scaled"]])) <= 1) {
    
    cat(sprintf("Fewer than 2 unique microbiome values. Diet: %s, Phenotype: %s, Feature: %s\n",
                this.diet, this.pheno, this.mb.feat))
    
  }
  stopifnot(length(unique(this.df$mb.scaled)) >= 2)
  
  # run models
  med.fit <- lmer(med.formula, data=this.df, REML=F)
  out.fit <- lmer(out.formula, data=this.df, REML=F)
  
  # do mediation
  this.med.res <- mediate(med.fit, out.fit,
                          treat="Diet.5mo.as.AL", treat.value=this.diet, control.value="AL",
                          mediator="mb.scaled", sims=1000)
  
  # save output as vector
  out.vec <- c(
    ACME_est=this.med.res$d0,
    ACME_2.5p=this.med.res$d0.ci[[1]],
    ACME_97.5p=this.med.res$d0.ci[[2]],
    ACME_p=this.med.res$d0.p,
    ADE_est=this.med.res$z0,
    ADE_2.5p=this.med.res$z0.ci[[1]],
    ADE_97.5p=this.med.res$z0.ci[[2]],
    ADE_p=this.med.res$z0.p,
    prop.med=this.med.res$n0,
    prop.med_2.5p=this.med.res$n0.ci[[1]],
    prop.med_97.5p=this.med.res$n0.ci[[2]],
    prop.med_p=this.med.res$n0.p,
    n=this.med.res$nobs,
    pheno=this.pheno,
    feature=this.mb.feat,
    diet=this.diet)
  
  return(out.vec)

  }

# Initialize cluster
n.cores <- n.threads
cat("Will use", n.cores, "cores...\n")
cl <- makeCluster(n.cores, outfile=log.path)
clusterExport(cl, c("run_mediation_one_pheno_one_mb_feat"))
registerDoParallel(cl)

# Import microbiome data
cat("Reading microbiome data...\n")
mb.df.feats.in.rows <- read.table(mb.data.path, sep="\t", header=T, row.names=1)

# Enumerate all features
all.mb.feats <- rownames(mb.df.feats.in.rows)

# Transpose to get stool.IDs in the rows
mb.df <- mb.df.feats.in.rows %>% t() %>% data.frame()

# Import phenotype data
cat("Reading phenotype data...\n")
pheno.df <- read.table(pheno.data.path, sep="\t", header=T)
phenos <- sort(colnames(pheno.df)[2:ncol(pheno.df)])

# Exclude SurvDays because it won't work with this longitudinal model
phenos <- str_subset(phenos, "SurvDays", negate=T)

# Import metadata
stool.meta.df <- read.table(stool.meta.path, sep="\t", header=T)
mouse.meta.df <- read.csv(mouse.meta.path)
stool.meta.annot.df <- stool.meta.df %>% 
  merge(mouse.meta.df, by.x="mouse.ID", by.y="MouseID")

# Identify pre-randomization samples as AL
stool.meta.annot.df <- stool.meta.annot.df %>%
  mutate(Diet.5mo.as.AL=factor(case_when(
    age.wks < 25 ~ "AL",
    TRUE ~ as.character(Diet)),
    levels=c("AL", "1D", "2D", "20", "40")))

# Combine microbiome, phenotypes, and metadata into one df
full.df <- pheno.df %>% 
  merge(stool.meta.annot.df, by="stool.ID") %>% 
  merge(mb.df %>% rownames_to_column("stool.ID"), by="stool.ID")

# Subset to this diet
this.diet.df <- full.df %>% 
  dplyr::filter(Diet.5mo.as.AL %in% c("AL", this.diet)) %>%
  mutate(Diet.5mo.as.AL=factor(Diet.5mo.as.AL, levels=c("AL", this.diet))) # convert to factor

# Subset to this phenotype
this.diet.this.pheno.df <- this.diet.df %>% 
  
  # remove missing phenotype values
  dplyr::filter(!is.na(.data[[this.pheno]])) %>%
  
  # scale phenotype and age
  mutate(pheno.scaled = as.numeric(scale(.data[[this.pheno]]))) %>%
  mutate(age.wks.scaled = as.numeric(scale(age.wks)))

# Loop over microbiome features
cat("Looping over microbiome features...\n")
this.diet.pheno.out.res <- foreach(
  this.mb.feat = all.mb.feats,
  .combine = "rbind",
  .packages = c("lme4", "mediation")) %dopar% {
    
    cat(this.mb.feat, "\n")

    # Run model for one microbiome feature
    this.mb.feat.out.list <- run_mediation_one_pheno_one_mb_feat(
      this.pheno = this.pheno,
      this.mb.feat = this.mb.feat,
      this.df = this.diet.this.pheno.df,
      this.diet = this.diet,
      med.formula = med.formula,
      out.formula = out.formula)
    
  } # microbiome features

# Write output for this diet-pheno  
out.df <- data.frame(this.diet.pheno.out.res)
saveRDS(out.df, out.path)
cat(sprintf("Wrote %s\n", out.path))
