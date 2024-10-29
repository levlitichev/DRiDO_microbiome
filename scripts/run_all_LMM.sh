# all ASReml microbiome models
sh run_all_asreml.sh # slow, cleans up intermediate files

# MaAsLin models
Rscript run_maaslin.R # does not clean up intermediate files

# we use lme4qtl just once for an Extended Data Figure, but it's very slow
Rscript run_lme4qtl.R # slow
Rscript collate_lme4qtl.R # does not clean up intermediate files

# microbiome-phenotype associations
Rscript run_lme4_mb_pheno_assoc.R
Rscript collate_lme4_mb_pheno_assoc.R

# cross-sectional microbiome-phenotype associations
# need to run the next script manually since it's .Rmd and not .R
# lm_cross_sectional_mb_pheno_assoc.Rmd

# N.B. mediation LMMs are run through Snakefile_mediation