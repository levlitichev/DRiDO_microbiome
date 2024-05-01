# Run a linear model for all microbiome-phenotype pairs with the following specification:

# y_pheno (scaled) ~ y_mb (f, scaled) + age (f, scaled) + diet (f) + mouse (r)

# Took ~10 min for ~100 microbiome features and 200 phenotypes.

# Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(lme4))
suppressPackageStartupMessages(library(foreach)) # for parallel for-loop
suppressPackageStartupMessages(library(doParallel)) # for parallel for-loop

### --- INPUTS --- ###

setwd("~/DRiDO_microbiome_github/scripts/")

out.dir <- "../results/lme4_mb_pheno_assoc_kraken/"
stopifnot(dir.exists(out.dir))

# Import microbiome data
cat("Reading microbiome data...\n")

# Genera
mb.df <- read.table(
  "../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt",
  sep="\t", header=T, row.names=1)
# mb.df[1:5, 1:3]

# Pathways
# mb.df.feats.in.rows <- read.table(
#   "../results/DO_pathway_log2tpm_filt_w_comm_n273x2997.txt",
#   sep="\t", header=T, row.names=1)
# mb.df <- data.frame(t(mb.df.feats.in.rows)) # transpose to get stool.ID in the rows
# mb.df[1:5, 1:3]

# Import phenotype data
cat("Reading phenotype data...\n")
pheno.df <- read.table(
  "../results/phenotype_values_closest_to_each_stool_ID_within_100_days_n2915x210_230811.txt", 
  sep="\t", header=T)

# Import metadata
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

# Define function that will work on one microbiome feature and one phenotype
run_lme4_one_pheno_one_mb_feat <- function(this.pheno, this.mb.feat, this.df, 
                                           full.formula, reduced.model) {
  
  # scale microbiome feature
  this.df[["mb.scaled"]] <- as.numeric(scale(this.df[[this.mb.feat]]))
  
  # full model
  full.model <- lmer(full.formula, data=this.df, REML=F)
  
  # LRT
  this.anova.res <- anova(full.model, reduced.model)
  
  # extract desired outputs
  out.list <- list(
    pheno = this.pheno,
    feature = this.mb.feat,
    mb.coef = fixef(full.model)[["mb.scaled"]],
    age.coef = fixef(full.model)[["age.wks.scaled"]],
    DR.1D.coef = fixef(full.model)[["Diet.5mo.as.AL1D"]],
    DR.2D.coef = fixef(full.model)[["Diet.5mo.as.AL2D"]],
    DR.20.coef = fixef(full.model)[["Diet.5mo.as.AL20"]],
    DR.40.coef = fixef(full.model)[["Diet.5mo.as.AL40"]],
    mb.pval = this.anova.res$`Pr(>Chisq)`[2]
  )
  
  return(out.list)
  
}

# Initialize cluster
foreach.out.log.path <- sprintf("%srun_lme4_foreach.log", out.dir)
cl <- makeCluster(6, outfile=foreach.out.log.path)
clusterExport(cl, c("run_lme4_one_pheno_one_mb_feat"))
registerDoParallel(cl)

# Keep track of features
all.mb.feats <- colnames(mb.df)
# length(all.mb.feats)

# Keep track of phenos
# Exclude SurvDays because it won't work with this longitudinal model
phenos <- sort(colnames(pheno.df)[2:ncol(pheno.df)])
phenos <- str_subset(phenos, "SurvDays", negate=T)
# length(phenos)

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
# dim(full.df)

# Define formulas for the linear models
reduced.formula <- as.formula("pheno.scaled ~ age.wks.scaled + Diet.5mo.as.AL + (1|mouse.ID)")
full.formula <- as.formula("pheno.scaled ~ mb.scaled + age.wks.scaled + Diet.5mo.as.AL + (1|mouse.ID)")

# Run all
cat("Looping over phenotypes...\n")
for (this.pheno in phenos) {
  
  cat(this.pheno, "\n")

  out.path <- sprintf("%s%s.Rds", out.dir, this.pheno)
  
  # Subset to this pheno
  this.df <- full.df %>% 
    
    # Remove missing phenotype values
    dplyr::filter(!is.na(.data[[this.pheno]]))

  # Make sure we have samples from each diet
  stopifnot(n_distinct(this.df$Diet.5mo.as.AL) == 5)
  
  # Scale phenotype and age
  this.df[["pheno.scaled"]] <- as.numeric(scale(this.df[[this.pheno]]))
  this.df[["age.wks.scaled"]] <- as.numeric(scale(this.df[["age.wks"]]))
  
  # Run model without microbiome (this step can be run just once per pheno)
  reduced.model <- lmer(reduced.formula, data=this.df, REML=F)

  # Loop over microbiome features
  cat("Looping over microbiome features...\n")
  this.pheno.out.res <- foreach(this.mb.feat = all.mb.feats,
			       .combine = "rbind",
			       .errorhandling = "remove",
			       .packages = c("lme4", "tidyverse")) %dopar% {

    # Run model for one microbiome feature
    this.mb.feat.out.list <- run_lme4_one_pheno_one_mb_feat(
      this.pheno = this.pheno,
      this.mb.feat = this.mb.feat,
      this.df = this.df,
      full.formula = full.formula,
      reduced.model = reduced.model)
    
  } # microbiome feature loop

  # Write output for this pheno  
  saveRDS(this.pheno.out.res, out.path)
  cat(sprintf("Wrote %s\n", out.path))

} # pheno loop