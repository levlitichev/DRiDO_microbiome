library(foreach)

in.paths <- list.files("output/", pattern="*.Rds", recursive=T, full.names=T)
out.path <- gzfile("mediation_n107_genus_feats_n209_phenos.txt.gz")

cat("Will combine", length(in.paths), "files\n")

# read in end combine all results
out.df <- foreach(this.path=in.paths, .combine="rbind") %do% {readRDS(this.path)}

# write output
write.table(out.df, out.path, sep="\t", row.names=F, quote=F)
