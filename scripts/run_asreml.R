# 2024-06-27: Script now accepts fixed and random formulas as arguments.
# Run a linear mixed model for all features.
# This script is fast enough to run locally, i.e. no need for a cluster (<5 minutes for 100 features and 3000 samples)

# load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(lubridate))
suppressPackageStartupMessages(library(asreml))
suppressPackageStartupMessages(library(foreach)) # for parallel for-loop
suppressPackageStartupMessages(library(doParallel)) # for parallel for-loop

# activating ASReml in the parallel for loop is a problem
asreml.license.offline(10)

# check available cores with detectCores()
cl <- makeCluster(6)
registerDoParallel(cl)

### --- COMMAND LINE INPUTS --- ###

# read in command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  for (ii in seq(length(args))) {cat(sprintf("args[%i]: %s\n", ii, args[ii]))}
  stop(sprintf(
    "Must provide 5 command line arguments. length(args) = %i",
    length(args)), call.=F)
}

# first argument is microbiome data, e.g. ../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt
mb.df <- read.table(args[1], sep="\t", header=T, row.names=1)
# mb.df[1:5, 1:5]

# second argument is output directory, e.g. ../results/asreml_kraken_genus/
# it must already exist
out.dir <- args[2]
if (!dir.exists(out.dir)) {
  stop(sprintf("Output directory must already exist. out.dir: %s", out.dir))
}

# third argument is the right hand side of the fixed effect formula
# e.g. "~ age.wks.scaled + Diet.5mo.as.AL + time.scaled"
fixef.formula.RHS <- args[3]
cat(sprintf("Fixed effects: %s\n", fixef.formula.RHS))

# fourth argument is the random effect formula
# e.g. "~ vm(Mouse, kinship.mat.x2) + ide(Mouse) + Cohort + Cage + Batch"
ranef.formula.str <- args[4]
cat(sprintf("Random effects: %s\n", ranef.formula.str))

# fifth argument is the random effect formula excluding genetics (this is needed to perform LRT)
# e.g. "~ ide(Mouse, kinship.mat.x2) + Cohort + Cage + Batch"
ranef.formula.no.genetics.str <- args[5]
cat(sprintf("Random effects without genetics: %s\n", ranef.formula.no.genetics.str))

### --- OTHER INPUTS --- ###

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

# create new columns related to date of stool collection
stool.meta.annot.df <- stool.meta.annot.df %>% 
  mutate(collection.date = ymd(date.stool.collection.approx)) %>% 
  mutate(quarter.collection.date = round_date(collection.date, "quarter")) %>% 
  mutate(days.from.first.stool.collection = time_length(
    collection.date - min(collection.date), unit="day"))

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
cat(sprintf("Number of samples after filtering: %i\n", nrow(final.mb.annot.df)))

# minor metadata tweaks
final.mb.annot.df <- final.mb.annot.df %>% 
  
  # pre-randomization timepoints should be considered AL
  mutate(Diet.5mo.as.AL = case_when(
    age.approx.months == 5 ~ "AL",
    TRUE ~ as.character(Diet))) %>% 
  
  # scale age and chronological time
  mutate(age.wks.scaled = scale(age.wks),
         time.scaled = scale(days.from.first.stool.collection)) %>% 
  
  # convert to factors
  mutate(Diet.5mo.as.AL=factor(Diet.5mo.as.AL, levels=c("AL", "1D", "2D", "20", "40")),
         Diet=factor(Diet, levels=c("AL", "1D", "2D", "20", "40")),
         Cohort=factor(Cohort),
         Cage=factor(Cage),
         Batch=factor(ext.batch),
         Mouse=factor(mouse.ID, levels=mice.to.keep),
         Quarter.Date=factor(quarter.collection.date))

# multiply kinship by 2
kinship.mat.x2 <- kinship.mat*2
rm(kinship.mat) # to make sure I don't use the wrong one

# define function to run on one feature at a time
run_one_feat_age_and_DR_fixed <- function(this.feat, this.df) {
  
  cat(this.feat, "\n")
  
  fixef.formula <- as.formula(paste(this.feat, fixef.formula.RHS))
  
  ranef.formula <- as.formula(ranef.formula.str)
  
  # fit model
  this.model <- asreml(
    fixed = fixef.formula,
    random = ranef.formula,
    data = this.df)
  
  # estimate standard error for genetic heritability
  # for vpredict, we need to create a formula that looks like "h2 ~ V2 / (V1 + V2 + V3)"
  # where V2 corresponds to the variance explaned by the genetics term
  num.ranef <- length(this.model$vparameters)
  stopifnot("vm(Mouse, kinship.mat.x2)" %in% names(this.model$vparameters))
  ii.for.genetics.ranef <- which(names(this.model$vparameters) == "vm(Mouse, kinship.mat.x2)")
  denom.for.vpredict <- paste0("V", seq(num.ranef), collapse=" + ")
  formula.str.for.vpredict <- sprintf("h2 ~ V%i / (%s)", ii.for.genetics.ranef, denom.for.vpredict)
  genetic.herit.se <- vpredict(this.model, as.formula(formula.str.for.vpredict))[["SE"]]
  
  # get fixed effect coefficients
  this.fixef.df <- data.frame(summary(this.model, coef=T)$coef.fixed)
  
  # run conditional Wald test (actually includes incremental result too)
  this.wald.df <- as.data.frame(wald(this.model, ssType = "conditional")$Wald)
  
  # run model without genetics
  # need to provide kinship matrix to ide even though we're not going to use it
  ranef.formula.no.genetics <- as.formula(paste(ranef.formula.no.genetics.str))
  
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

