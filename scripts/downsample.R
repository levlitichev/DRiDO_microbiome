setwd("~/DRiDO_microbiome_github/scripts/")

# import genus-level data
genus.data.df <- read.table(
  "../results/kraken_genus_clr_filt_w_comm_n107x2997_231017.txt", 
  sep="\t", header=T, row.names=1)

# get 519 random stool IDs so that n is equal to per-age analysis
random.n.stool.IDs <- scan(
  "../results/stool_IDs_for_asreml_downsampled_n519.txt",
  what="character")
length(random.n.stool.IDs)

# downsample
genus.data.df.downsampled <- genus.data.df[, random.n.stool.IDs]

# export
# genus.data.df.downsampled %>% 
#   rownames_to_column("feature") %>%
#   write.table(sprintf(
#     "../results/kraken_genus_clr_filt_w_comm_downsampled_n%ix%i.txt",
#     nrow(genus.data.df.downsampled), ncol(genus.data.df.downsampled)),
#     sep="\t", quote=F, row.names=F)
