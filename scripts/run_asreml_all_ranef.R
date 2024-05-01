# Run a linear model for all features with the following specification (all random effects):

# y_mb ~ age (r) + diet (r) + [genetics] (r) + mouse (r) + batch (r) + cohort (r) + cage (r) 

# pretty slow, takes ~1 hour for ~100 features

# load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(asreml))
suppressPackageStartupMessages(library(foreach)) # for parallel for-loop
suppressPackageStartupMessages(library(doParallel)) # for parallel for-loop

# activating ASReml in the parallel for loop is a problem
asreml.license.offline(10)

# change working directory
setwd("~/DRiDO_microbiome_github/scripts/")

# check available cores with detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)

### --- INPUTS --- ###

# specify output directory
out.dir <- "../results/asreml_kraken_genus_all_ranef/"
stopifnot(dir.exists(out.dir))

# import microbiome data
mb.df <- read.table(
  "../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt",
  # "../results/DO_pathway_log2tpm_filt_w_comm_n273x2997.txt",
  sep="\t", header=T, row.names=1)
# mb.df[1:5, 1:3]

# import kinship matrix
kinship.df <- read.csv(
  "../data/kinship.all_chroms_downloaded_from_wright22.csv",
  check.names=F) %>% 
  column_to_rownames("106953") # don't know why the first column is called this...

# import metadata
stool.meta.df <- read.table(
  "../data/metadata/stool_metadata_after_QC_no_controls_n2997_240418.txt", 
  sep="\t", header=T)
mouse.meta.df <- read.csv(
  "../data/metadata/AnimalData_Processed_20230712.csv") %>% 
  mutate(Cage = paste0("c", HID))

stool.meta.annot.df <- merge(
  stool.meta.df, mouse.meta.df, 
  by.x="mouse.ID", by.y="MouseID")

### --- END INPUTS --- ###

# keep track of features
all.feats <- rownames(mb.df)
cat(sprintf("Running ASReml on %i features...\n", length(all.feats)))

# scale features
mb.scaled.df <- mb.df %>% 
  t() %>% # scale expects samples in the rows
  scale() %>% 
  data.frame() %>% # keep samples in the rows
  rownames_to_column("stool.ID")
# dim(mb.scaled.df)

# add metadata to microbiome data
mb.scaled.annot.df <- merge(
  mb.scaled.df,
  stool.meta.annot.df,
  by="stool.ID")
# dim(mb.scaled.annot.df)

# subset to mice with both microbiome and genotype information
mice.to.keep <- intersect(
  unique(mb.scaled.annot.df$mouse.ID),
  rownames(kinship.df)
)
# length(mice.to.keep)

# subset kinship and convert to matrix
kinship.mat <- as.matrix(kinship.df[mice.to.keep, mice.to.keep])
# dim(kinship.mat)

# subset microbiome data
final.mb.annot.df <- mb.scaled.annot.df %>% 
  dplyr::filter(mouse.ID %in% mice.to.keep)
# dim(final.mb.annot.df)

# minor metadata tweaks
final.mb.annot.df <- final.mb.annot.df %>% 
  
  # pre-randomization timepoints should be considered AL
  mutate(Diet.5mo.as.AL = case_when(
    age.approx.months == 5 ~ "AL",
    TRUE ~ as.character(Diet))) %>% 
  
  # convert age to factor to use as random effect
  mutate(Age=factor(age.approx.months)) %>% 
  
  # convert to factors
  mutate(Diet.5mo.as.AL=factor(Diet.5mo.as.AL, levels=c("AL", "1D", "2D", "20", "40")),
         Diet=factor(Diet, levels=c("AL", "1D", "2D", "20", "40")),
         Cohort=factor(Cohort),
         Cage=factor(Cage),
         Batch=factor(ext.batch),
         Mouse=factor(mouse.ID, levels=mice.to.keep))

# multiply kinship by 2
kinship.mat.x2 <- kinship.mat*2
rm(kinship.mat) # to make sure I don't use the wrong one

# define function to run on one feature at a time
run_one_feat_age_and_DR_ranef <- function(this.feat, this.df) {
  
  cat(this.feat, "\n")
  
  fixef.formula <- as.formula(paste(
    this.feat, "~ 1"
  ))
  
  ranef.formula <- as.formula(paste(
    "~ Age + Diet.5mo.as.AL + vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch"
  ))
  
  # ~~~ FULL MODEL ~~~
  this.model <- asreml(
    fixed = fixef.formula,
    random = ranef.formula,
    data = this.df)
  
  # ~~~ NO AGE ~~~
  ranef.formula.no.age <- as.formula(paste(
    "~ Diet.5mo.as.AL + vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch"
  ))
  
  this.model.no.age <- asreml(
    fixed = fixef.formula,
    random = ranef.formula.no.age,
    data = this.df)
  
  this.LRT.pval.age <- 1 - pchisq(2 * (this.model$loglik - this.model.no.age$loglik), 1)
  
  # ~~~ NO DIET ~~~
  ranef.formula.no.diet <- as.formula(paste(
    "~ Age + vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch"
  ))
  
  this.model.no.diet <- asreml(
    fixed = fixef.formula,
    random = ranef.formula.no.diet,
    data = this.df)
  
  this.LRT.pval.diet <- 1 - pchisq(2 * (this.model$loglik - this.model.no.diet$loglik), 1)
  
  # ~~~ NO MOUSE ~~~
  ranef.formula.no.mouse <- as.formula(paste(
    "~ Age + Diet.5mo.as.AL + vm(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
  ))
  
  this.model.no.mouse <- asreml(
    fixed = fixef.formula,
    random = ranef.formula.no.mouse,
    data = this.df)
  
  this.LRT.pval.mouse <- 1 - pchisq(2 * (this.model$loglik - this.model.no.mouse$loglik), 1)
  
  # ~~~ NO GENETICS ~~~
  ranef.formula.no.genetics <- as.formula(paste(
    "~ Age + Diet.5mo.as.AL + ide(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
  ))
  
  this.model.no.genetics <- asreml(
    fixed = fixef.formula,
    random = ranef.formula.no.genetics,
    data = this.df)
  
  this.LRT.pval.genetics <- 1 - pchisq(2 * (this.model$loglik - this.model.no.genetics$loglik), 1)
  
  # ~~~ NO COHORT ~~~
  ranef.formula.no.cohort <- as.formula(paste(
    "~ Age + Diet.5mo.as.AL + vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cage + Batch"
  ))
  
  this.model.no.cohort <- asreml(
    fixed = fixef.formula,
    random = ranef.formula.no.cohort,
    data = this.df)
  
  this.LRT.pval.cohort <- 1 - pchisq(2 * (this.model$loglik - this.model.no.cohort$loglik), 1)
  
  # ~~~ NO CAGE ~~~
  ranef.formula.no.cage <- as.formula(paste(
    "~ Age + Diet.5mo.as.AL + vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Batch"
  ))
  
  this.model.no.cage <- asreml(
    fixed = fixef.formula,
    random = ranef.formula.no.cage,
    data = this.df)
  
  this.LRT.pval.cage <- 1 - pchisq(2 * (this.model$loglik - this.model.no.cage$loglik), 1)
  
  # ~~~ NO BATCH ~~~
  ranef.formula.no.batch <- as.formula(paste(
    "~ Age + Diet.5mo.as.AL + vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage"
  ))
  
  this.model.no.batch <- asreml(
    fixed = fixef.formula,
    random = ranef.formula.no.batch,
    data = this.df)
  
  this.LRT.pval.batch <- 1 - pchisq(2 * (this.model$loglik - this.model.no.batch$loglik), 1)
  
  # ~~~ CREATE OUTPUT DF ~~~
  
  # calculate PVE
  this.varcomp.df <- summary(this.model)$varcomp %>%
    mutate(PVE=component/(sum(component)))
  
  # add column for feature
  this.varcomp.df$feature <- this.feat
  
  # move rownames to new column
  this.varcomp.df <- this.varcomp.df %>% 
    rownames_to_column("term")
  
  # make cleaner term column
  this.varcomp.df <- this.varcomp.df %>% 
    mutate(term.clean=case_when(
      term == "Diet.5mo.as.AL" ~ "Diet",
      term == "vm(Mouse, kinship.mat.x2)" ~ "Genetics",
      term == "ide(Mouse)" ~ "Mouse",
      term == "units!R" ~ "Residual",
      TRUE ~ term))
  
  # convert p-values to dataframe
  LRT.pval.df <- data.frame(LRT.pval = c(
    Age = this.LRT.pval.age,
    Diet = this.LRT.pval.diet,
    Mouse = this.LRT.pval.mouse,
    Genetics = this.LRT.pval.genetics,
    Cohort = this.LRT.pval.cohort,
    Cage = this.LRT.pval.cage,
    Batch = this.LRT.pval.batch
  )) %>% rownames_to_column("term.clean")
  
  # merge p-values with random effect dataframe
  out.df <- this.varcomp.df %>% 
    merge(LRT.pval.df, by="term.clean", all.x=T)
  
  # write output to file
  out.path <- sprintf("%s%s.txt", out.dir, this.feat)
  cat(sprintf("Writing output to %s\n", out.path))
  write.table(out.df, file=out.path, sep="\t", row.names=F, quote=F)
  
}

# run all
asreml.options(trace=F, workspace="1gb")

IGNORE <- foreach(
  this.feat = all.feats,
  .packages = c("tidyverse", "asreml")) %dopar% {
    
    .GlobalEnv$kinship.mat.x2 <- kinship.mat.x2
    run_one_feat_age_and_DR_ranef(
      this.feat, 
      this.df = final.mb.annot.df)
  }

cat("Done!\n")