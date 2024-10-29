in.dir="~/DRiDO_microbiome_github/results/lme4qtl_kraken_genus_w_time_ranef/"
out.path="~/DRiDO_microbiome_github/results/lme4qtl_kraken_genus_w_time_ranef/herit.txt"

library(tidyverse)
library(foreach)

# vector of all lme4qtl result text files
lme4qtl.paths.vec <- list.files(in.dir, ".txt", full.names=T)

lme4qtl.herit.df <- foreach(this.path = lme4qtl.paths.vec, .combine="rbind") %do% {
  
  # import each lme4qtl result as df
  one.lme4qtl.res.df <- read.table(this.path, sep="\t", header=T)
  
  # subset to just the row that contains heritability
  one.lme4qtl.herit.df <- one.lme4qtl.res.df %>% 
    dplyr::filter(grp == "mouse.ID")
  
  return(one.lme4qtl.herit.df)
}

# write output
lme4qtl.herit.df %>% 
  write.table(out.path, sep="\t", row.names=F, quote=F)
