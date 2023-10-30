library(tidyverse)
library(foreach)

list.of.RDS.paths <- list.files("../results/lme4_mb_pheno_assoc_kraken/", pattern="*.Rds", full.names=T)

# check if all runs finished
completed.phenos <- str_replace(basename(list.of.RDS.paths), ".Rds", "")
length(completed.phenos) # 209 is correct, because SurvDays won't work

# loop over all Rds files
out.df <- foreach(this.path = list.of.RDS.paths, .combine="rbind") %do% {
  
  this.res <- readRDS(this.path)
  
  data.frame(this.res) %>% 
    unnest(everything())
  
}

# summarise number of microbiome features and phenotypes
n.mb.feats <- n_distinct(out.df$feature)
n.phenos <- n_distinct(out.df$pheno)

# export
out.df %>% 
  write.table(
    gzfile(sprintf("../results/lme4_mb_pheno_assoc_n%i_kraken_feats_n%i_phenos.txt.gz", n.mb.feats, n.phenos)),
    sep="\t", quote=F, row.names=F)
