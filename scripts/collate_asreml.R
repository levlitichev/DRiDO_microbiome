suppressPackageStartupMessages(library(tidyverse))

# read command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop(sprintf(
    "Must provide only one command line argument: a directory with RDS files. args: %s",
    args), call.=F)
}

# the only input should be the output of run_asreml*.R: a directory with RDS files, e.g. ../results/asreml_kraken_genus/
dir.w.RDS.files <- args[1]
cat(sprintf("Input directory: %s\n", dir.w.RDS.files))

# get list of RDS files to collate
list.of.RDS.paths <- list.files(dir.w.RDS.files, pattern=".Rds", full.names=T)
cat(sprintf("Collating %i results...\n", length(list.of.RDS.paths)))

# initialize outputs
list.of.fixed.dfs <- vector("list", length(list.of.RDS.paths))
list.of.wald.dfs <- vector("list", length(list.of.RDS.paths))
list.of.varcomp.dfs <- vector("list", length(list.of.RDS.paths))
feat.name.vec <- rep("", length(list.of.RDS.paths))
pval.vec <- rep(0, length(list.of.RDS.paths))
herit.se.vec <- rep(0, length(list.of.RDS.paths))

# loop over RDS files
for (ii in seq(length(list.of.RDS.paths))) {
  
  this.res <- readRDS(list.of.RDS.paths[[ii]])
  this.feat <- tools::file_path_sans_ext(basename(list.of.RDS.paths[[ii]]))
  
  list.of.fixed.dfs[[ii]] <- this.res[["fixed.df"]]
  list.of.wald.dfs[[ii]] <- this.res[["wald.df"]]
  list.of.varcomp.dfs[[ii]] <- this.res[["varcomp.df"]]
  pval.vec[[ii]] <- this.res[["genetics.LRT.pval"]]
  herit.se.vec[[ii]] <- this.res[["genetic.herit.se"]]
  feat.name.vec[[ii]] <- this.feat
  
}

# concatenate
fixed.df <- do.call(rbind, list.of.fixed.dfs)
wald.df <- do.call(rbind, list.of.wald.dfs)
varcomp.df <- do.call(rbind, list.of.varcomp.dfs)
pval.df <- data.frame(
  feature = feat.name.vec,
  genetics.LRT.pval = pval.vec)
herit.se.df <- data.frame(
  feature = feat.name.vec,
  herit.se = herit.se.vec)

# combine heritability estimate, heritability standard error, and LRT p-values
herit.df <- merge(
  varcomp.df %>% 
    dplyr::filter(term == "vm(Mouse, kinship.mat.x2)") %>% 
    dplyr::select(feature, PVE),
  herit.se.df, by="feature") %>% 
  merge(pval.df, by="feature") %>% 
  dplyr::rename(herit = PVE,
                LRT.pval = genetics.LRT.pval)

# combine fixed effect coefficients in fixed.df with Wald test p-values
# N.B. wald.df has one p-value for the diet coefficient, but fixed.df has a separate coefficient for each diet
out.fixef.df <- merge(
  
  fixed.df %>% 
    separate(col=term, into=c("fixef", "fixef.value"),  sep="_", fill="right"),
  
  wald.df %>% 
    dplyr::rename(fixef = term),
  
  by=c("fixef", "feature"))

# confirm the sizes are what we expect  
# dim(wald.df)
# dim(fixed.df)
# dim(out.fixef.df)

# export fixef.df, varcomp.df (AKA ranef.df), and herit.df
write.table(out.fixef.df, sprintf("%s/fixef.txt", dir.w.RDS.files),
            sep="\t", quote=F, row.names=F)

write.table(varcomp.df, sprintf("%s/ranef.txt", dir.w.RDS.files),
            sep="\t", quote=F, row.names=F)

write.table(herit.df, sprintf("%s/herit.txt", dir.w.RDS.files),
            sep="\t", quote=F, row.names=F)

cat("Wrote fixef.txt, ranef.txt, and herit.txt.\n")