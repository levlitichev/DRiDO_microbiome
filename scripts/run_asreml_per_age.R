# Run a linear model for all features with the following specification:

# After randomization:  y_mb ~ diet (f) + [genetics] (r) + batch (r) + cohort (r) + cage (r) 
# At 5 months:          y_mb ~ [genetics] (r) + batch (r) + cohort (r) + cage (r) 

# this script runs quickly, just a minute or two per age on 100 features

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
out.dir <- "../results/asreml_kraken_genus_per_age/"
stopifnot(dir.exists(out.dir))

# 5, 10, 16, 22, 28
# remember to also change the formula in the function below
this.age <- 28
cat(this.age, "months\n")

# import microbiome data
mb.df <- read.table(
  "../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt",
  sep="\t", header=T, row.names=1)
# mb.df[1:5, 1:3]

# import kinship matrix
kinship.df <- read.csv(
  "../data/kinship.all_chroms_downloaded_from_wright22.csv",
  check.names=F) %>% 
  column_to_rownames("106953") # don't know why the first column is called this...

# import metadata
stool.meta.df <- read.table(
  "../data/metadata/stool_metadata_after_QC_no_controls_n2997_230620.txt", 
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

# get mice with microbiome at this age
# (sort to help make sure kinship matrix and microbiome data match up later on)
this.age.mice <- stool.meta.annot.df %>% 
  dplyr::filter(age.approx.months == this.age) %>% 
  arrange(mouse.ID) %>% 
  pull(mouse.ID)
# length(this.age.mice)

# subset to mice with both microbiome at this age and genotype information
this.age.mice.to.keep <- intersect(
  this.age.mice,
  rownames(kinship.df)
)
# length(this.age.mice.to.keep)

# subset metadata to this age and these mice
this.age.stool.meta.annot.df <- stool.meta.annot.df %>% 
  dplyr::filter(age.approx.months == this.age) %>% 
  dplyr::filter(mouse.ID %in% this.age.mice.to.keep)
# nrow(this.age.stool.meta.annot.df)

# subset microbiome data to these mice and this age
this.age.mb.df <- mb.df[, this.age.stool.meta.annot.df$stool.ID]
# dim(this.age.mb.df)

# scale features
this.age.mb.scaled.df <- this.age.mb.df %>% 
  t() %>% # transpose to get samples in the rows
  scale() %>% 
  data.frame() %>% # keep samples in the rows
  rownames_to_column("stool.ID")
# dim(this.age.mb.scaled.df)

# add metadata to microbiome data
this.age.mb.scaled.annot.df <- merge(
  this.age.mb.scaled.df,
  this.age.stool.meta.annot.df,
  by="stool.ID")
# dim(this.age.mb.scaled.annot.df)

# make sure mice in mb df are in same order as mice in kinship matrix
stopifnot(identical(this.age.mb.scaled.annot.df$mouse.ID, this.age.mice.to.keep))

# subset kinship and convert to matrix
this.age.kinship.mat <- as.matrix(kinship.df[this.age.mice.to.keep, this.age.mice.to.keep])
# dim(this.age.kinship.mat)

# minor metadata tweaks
final.mb.annot.df <- this.age.mb.scaled.annot.df %>% 
  
  # pre-randomization timepoints should be considered AL
  mutate(Diet.5mo.as.AL = case_when(
    age.approx.months == 5 ~ "AL",
    TRUE ~ as.character(Diet))) %>% 
  
  # convert to factors
  mutate(Diet.5mo.as.AL=factor(Diet.5mo.as.AL, levels=c("AL", "1D", "2D", "20", "40")),
         Diet=factor(Diet, levels=c("AL", "1D", "2D", "20", "40")),
         Cohort=factor(Cohort),
         Cage=factor(Cage),
         Batch=factor(ext.batch),
         Mouse=factor(mouse.ID, levels=this.age.mice.to.keep))

# multiply kinship by 2
kinship.mat.x2 <- this.age.kinship.mat*2
rm(this.age.kinship.mat) # to make sure I don't use the wrong one

# only consider microbiome features with enough unique values?
# final.mb.annot.df[, all.feats] %>% 
#   apply(MARGIN=2, FUN=n_distinct) %>% 
#   sort() %>% head()

# define function to run on one feature at a time
run_one_feat_one_age_DR_fixed <- function(this.feat, this.df) {
  
  cat(this.feat, "\n")
  
  fixef.formula <- as.formula(paste(
    # this.feat, "~ 1" # 5 months
    this.feat, "~ Diet.5mo.as.AL" # after start of DR
  ))
  
  ranef.formula <- as.formula(paste(
    "~ vm(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
  ))
  
  # fit model
  this.model <- asreml(
    fixed = fixef.formula,
    random = ranef.formula,
    data = this.df)
  
  # estimate standard error for genetic heritability
  genetic.herit.se <- vpredict(this.model, h2 ~ V4 / (V1 + V2 + V3 + V4 + V5))[["SE"]]
  
  # get fixed effect coefficients
  this.fixef.df <- data.frame(summary(this.model, coef=T)$coef.fixed)
  
  # run conditional Wald test (actually includes incremental result too)
  this.wald.df <- as.data.frame(wald(this.model, ssType = "conditional")$Wald)
  
  # run model without genetics
  # need to provide kinship matrix to ide even though we're not going to use it
  ranef.formula.no.genetics <- as.formula(paste(
    "~ Cohort + Cage + Batch"
  ))
  
  this.model.no.genetics <- asreml(
    fixed = fixef.formula,
    random = ranef.formula.no.genetics,
    data = this.df)
  
  # LRT p-value for genetics
  this.genetics.LRT.pval <- 1 - pchisq(2 * (this.model$loglik - this.model.no.genetics$loglik), 1)
  
  # output random effect dataframe
  this.varcomp.df <- summary(this.model)$varcomp %>%
    
    # calculate PVE
    mutate(PVE=component/(sum(component)))
  
  # add column for feature to all dataframes
  this.varcomp.df$feature <- this.feat
  this.fixef.df$feature <- this.feat
  this.wald.df$feature <- this.feat
  
  # Wald df has two columns called `Pr(chisq)`; need to fix this for rownames_to_column to work
  colnames(this.wald.df)[3] <- "Pr(chisq)(inc)"
  colnames(this.wald.df)[5] <- "Pr(chisq)(con)"
  
  # save output as Rds
  out.list <- list(
    fixed.df=this.fixef.df %>% rownames_to_column("term"),
    wald.df=this.wald.df %>% rownames_to_column("term"),
    varcomp.df=this.varcomp.df %>% rownames_to_column("term"),
    genetics.LRT.pval=this.genetics.LRT.pval,
    genetic.herit.se=genetic.herit.se
    )
  
  # write output to file
  out.path <- sprintf("%sM%i_%s.Rds", out.dir, this.age, this.feat)
  cat(sprintf("Writing output to %s\n", out.path))
  saveRDS(out.list, file=out.path)
  
}

# run all
asreml.options(trace=F, workspace="1gb")

IGNORE <- foreach(
  this.feat = all.feats, 
  .packages = c("tidyverse", "asreml")) %dopar% {
    
    .GlobalEnv$kinship.mat.x2 <- kinship.mat.x2
    run_one_feat_one_age_DR_fixed(
      this.feat, 
      this.df=final.mb.annot.df)
  }

cat("Done!\n")

