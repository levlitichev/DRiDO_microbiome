# same as run_asreml.R except with downsampling

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
out.dir <- "../results/asreml_kraken_genus_downsampled/"
stopifnot(dir.exists(out.dir))

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
  "../data/metadata/stool_metadata_after_QC_no_controls_n2997_240418.txt", 
  sep="\t", header=T)
mouse.meta.df <- read.csv(
  "../data/metadata/AnimalData_Processed_20230712.csv") %>% 
  mutate(Cage = paste0("c", HID))

stool.meta.annot.df <- merge(
  stool.meta.df, mouse.meta.df, 
  by.x="mouse.ID", by.y="MouseID")

# subset to ~519 samples so that n is equal to cross-sectional analysis
random.n.stool.IDs <- scan(
  "../results/stool_IDs_for_asreml_downsampled_n519.txt",
  what="character")
length(random.n.stool.IDs)

### --- END INPUTS --- ###

# downsample
mb.df <- mb.df[, random.n.stool.IDs]

# keep track of features
all.feats <- rownames(mb.df)
cat(sprintf("Running ASReml on %i features...\n", length(all.feats)))

# scale features
mb.scaled.df <- mb.df %>% 
  t() %>% # transpose to get samples in the rows
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
  
  # scale age
  mutate(age.wks.scaled=scale(age.wks)) %>% 
  
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
run_one_feat_age_and_DR_fixed <- function(this.feat, this.df) {
  
  cat(this.feat, "\n")
  
  fixef.formula <- as.formula(paste(
    this.feat, "~ age.wks.scaled + Diet.5mo.as.AL"
  ))
  
  ranef.formula <- as.formula(paste(
    "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch"
  ))
  
  # fit model
  this.model <- asreml(
    fixed = fixef.formula,
    random = ranef.formula,
    data = this.df)
  
  # estimate standard error for genetic heritability
  genetic.herit.se <- vpredict(this.model, h2 ~ V4 / (V1 + V2 + V3 + V4 + V5 + V6))[["SE"]]
  
  # get fixed effect coefficients
  this.fixef.df <- data.frame(summary(this.model, coef=T)$coef.fixed)
  
  # run conditional Wald test (actually includes incremental result too)
  this.wald.df <- as.data.frame(wald(this.model, ssType = "conditional")$Wald)
  
  # run model without genetics
  # need to provide kinship matrix to ide even though we're not going to use it
  ranef.formula.no.genetics <- as.formula(paste(
    "~ ide(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
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
    genetic.herit.se=genetic.herit.se,
    n=nrow(this.df)
    )
  
  # write output to file
  out.path <- sprintf("%s%s.Rds", out.dir, this.feat)
  cat(sprintf("Writing output to %s\n", out.path))
  saveRDS(out.list, file=out.path)
  
}

# run all
asreml.options(trace=F, workspace="1gb")

IGNORE <- foreach(
  this.feat = all.feats, 
  .packages = c("tidyverse", "asreml")) %dopar% {
    
    .GlobalEnv$kinship.mat.x2 <- kinship.mat.x2
    run_one_feat_age_and_DR_fixed(
      this.feat, 
      this.df=final.mb.annot.df)
  }

cat("Done!\n")

