# Run a linear mixed model for all features.
# Very similar to run_asreml.R but uses lme4qtl instead of ASReml.
# This script takes > 1 hour on a laptop (100 features and 3000 samples)

# set working directory
setwd("~/DRiDO_microbiome_github/scripts")

# load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(lme4qtl))
suppressPackageStartupMessages(library(foreach)) # for parallel for-loop
suppressPackageStartupMessages(library(doParallel)) # for parallel for-loop

### --- INPUTS --- ###

# import microbiome data
mb.df <- read.table(
  "../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt",
  sep="\t", header=T, row.names=1)
# mb.df[1:5, 1:5]

# make sure output directory already exists
out.dir <- "../results/lme4qtl_kraken_genus/"
if (!dir.exists(out.dir)) {
  stop(sprintf("Output directory must already exist. out.dir: %s", out.dir))
}

# specify model
# 1|mouse.ID accounts for additive genetic effects, while
# 1|mouse.ID.repeatability accounts for repeatability
formula.RHS <- "~ age.wks.scaled + Diet.5mo.as.AL + (1|Quarter.Date) + (1|mouse.ID) + (1|mouse.ID.repeatability) + (1|Batch) + (1|Cohort) + (1|Cage)"
cat(sprintf("Model: %s\n", formula.RHS))

# specify model without genetics to calculate LRT
formula.no.genetics.RHS <- "~ age.wks.scaled + Diet.5mo.as.AL + (1|Quarter.Date) + (1|mouse.ID.repeatability) + (1|Batch) + (1|Cohort) + (1|Cage)"
cat(sprintf("Model without genetics: %s\n", formula.no.genetics.RHS))

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

# define function to run on one feature at a time
run_lme4qtl_one_feat <- function(this.feat, this.df, this.kinship.mat) {
  
  cat(this.feat, "\n")

  # define model
  this.formula <- as.formula(paste(this.feat, formula.RHS))
  
  # fit model
  this.model <- relmatLmer(
    this.formula,
    this.df,
    relmat = list(mouse.ID = this.kinship.mat))
  
  # can get 95% confidence interval with the following code, but it's super slow
  # prof <- profile(this.lme4qtl.res)
  
  # get random effects
  ranef.df <- VarProp(this.model) %>% 
    mutate(feature=this.feat)
  
  # define model without genetics
  this.reduced.formula <- as.formula(paste(this.feat, formula.no.genetics.RHS))
  
  # I'm using relmatLmer, but lmer would also work
  this.model.no.genetics <- relmatLmer(
    this.reduced.formula,
    this.df)
  
  # LRT p-value for genetics
  full.loglik <- logLik(this.model, REML=T)[[1]]
  reduced.loglik <- logLik(this.model.no.genetics, REML=T)[[1]]
  this.genetics.LRT.pval <- 1 - pchisq(2 * (full.loglik - reduced.loglik), 1)
  
  # add genetics p-value 
  ranef.df$LRT.pval <- NA
  ranef.df[ranef.df$grp == "mouse.ID", "LRT.pval"] <- this.genetics.LRT.pval
  
  # write output to file
  out.path <- sprintf("%s%s.txt", out.dir, this.feat)
  cat(sprintf("Writing output to %s\n", out.path))
  write.table(ranef.df, file=out.path, sep="\t", quote=F, row.names=F)
  
}

# initialize cluster
foreach.out.log.path <- sprintf("%srun_lme4qtl_foreach.log", out.dir)
cl <- makeCluster(6, outfile=foreach.out.log.path)
clusterExport(cl, c("run_lme4qtl_one_feat"))
registerDoParallel(cl)

# keep track of features
all.feats <- rownames(mb.df)
cat(sprintf("Running lme4qtl on %i features...\n", length(all.feats)))

# scale features
mb.scaled.df <- mb.df %>% 
  t() %>% # scale expects samples in the rows
  scale() %>% 
  data.frame() %>% # keep samples in the rows
  rownames_to_column("stool.ID")
# dim(mb.scaled.df)

# create new columns related to date of stool collection
stool.meta.annot.df <- stool.meta.annot.df %>% 
  mutate(collection.date = ymd(date.stool.collection.approx)) %>% 
  mutate(quarter.collection.date = round_date(collection.date, "quarter"))

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
  
  # scale age
  mutate(age.wks.scaled=scale(age.wks)) %>% 
  
  # convert to factors
  mutate(Diet.5mo.as.AL=factor(Diet.5mo.as.AL, levels=c("AL", "1D", "2D", "20", "40")),
         Diet=factor(Diet, levels=c("AL", "1D", "2D", "20", "40")),
         Cohort=factor(Cohort),
         Cage=factor(Cage),
         Batch=factor(ext.batch),
         mouse.ID.repeatability=mouse.ID,
         Quarter.Date=factor(quarter.collection.date))

# multiply kinship by 2
kinship.mat.x2 <- kinship.mat*2
rm(kinship.mat) # to make sure I don't use the wrong one

# run all
IGNORE <- foreach(
  this.feat = all.feats, 
  .packages = c("tidyverse", "lme4qtl")) %dopar% {
    
    .GlobalEnv$kinship.mat.x2 <- kinship.mat.x2
    run_lme4qtl_one_feat(
      this.feat, 
      this.df=final.mb.annot.df,
      this.kinship.mat=kinship.mat.x2)
  }

cat("Done!\n")
