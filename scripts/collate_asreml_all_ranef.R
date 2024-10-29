suppressPackageStartupMessages(library(foreach))

# read command line arguments
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  stop(sprintf(
    "Must provide only one command line argument: a directory with .txt files. args: %s",
    args), call.=F)
}

# the only input should be the output of run_asreml_all_ranef*.R: a directory with .txt files, e.g. ../results/asreml_kraken_genus_all_ranef_w_time/
dir.w.txt.files <- args[1]
cat(sprintf("Input directory: %s\n", dir.w.txt.files))

# get list of text files to collate
list.of.txt.paths <- list.files(dir.w.txt.files, pattern=".txt", full.names=T)
cat(sprintf("Collating %i results...\n", length(list.of.txt.paths)))

# use rbind to collate
ranef.df.all.ranef.model <- foreach(this.path=list.of.txt.paths, .combine="rbind") %do% {
  read.table(this.path, sep="\t", header=T)
}

# write output
out.path <- file.path(dir.w.txt.files, "ranef.txt")
write.table(
  ranef.df.all.ranef.model,
  out.path, sep="\t", quote=F, row.names=F)
cat(sprintf("Wrote %s\n", out.path))