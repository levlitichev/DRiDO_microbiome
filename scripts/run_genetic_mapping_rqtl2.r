#!/usr/bin/env Rscript

suppressMessages(library(qtl2))
suppressMessages(library(qtl2convert))
suppressMessages(library(stringr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

####################################################################
###### set inputs and working directory

setwd("~/DRiDO_microbiome_github/")

geno_probability_file <- c('../data/prob.8state.allele.qtl2_200131.Rdata')     #obtain from FigShare XXXX

map_file <- c('../data/gm_uwisc_v1_filter_200129.pmap.csv')  #from https://github.com/kbroman/MUGAarrays/tree/main/UWisc

query_variants_file <- c('../data/cc_variants.sqlite') #from: https://figshare.com/articles/dataset/SQLite_database_of_variants_in_Collaborative_Cross_founder_mouse_strains/5280229   

query_genes_file <- c('../data/mouse_genes_mgi.sqlite')   #from: https://figshare.com/articles/dataset/SQLite_database_with_MGI_mouse_gene_annotations_from_Mouse_Genome_Informatics_MGI_at_The_Jackson_Laboratory/5286019

pheno_data_file <- c('../data/all_genus_features_clr_mice_in_rows_n910x535_230816.edit.txt')

meta_data_file  <- c('../data/jax_cr_genwave_191122.csv')

###### user to edit
focal_phenotype <- c('M10_Roseburia') #user to edit 
output_base <- paste("../results/rqtl2/",focal_phenotype,sep="")

####################################################################
## initialize varaibles 
##
n_reps <- 100 ## number of permutations to run to assess statistical significance of QTL peaks
model_type <- c('add')  ##type of genetic model to run - additive
focal_additive_cov <- c('cm_DietGenwave')  ##fixed effect covariates - diet + generation_wave

####################################################################
## LOAD.DATA_1 load genotype probability file
load(geno_probability_file)

# Edit name format from "Calico_Life_Sciences_Freund_MURGIGV01_20180621_DO.2D.4170_A1" to "DO-2D-4170".
print("rename mouse_ids")
mouse_IDs <- dimnames(apr$`1`)[[1]] %>% substr(48,57) %>% gsub(pattern = '[.]', replacement = '-')
if(nchar(mouse_IDs[1]) > 1) { 
  for(i in names(apr)){
    dimnames(apr[[i]])[[1]] <- mouse_IDs
  }   
}
print("edit marker names")
for(i in names(apr)){
  dimnames(apr[[i]])[[3]] <- dimnames(apr[[i]])[[3]] %>% gsub(pattern = '-', replacement = '.')
}

print('number of samples with genotype data')
print(length(mouse_IDs))

####################################################################
## LOAD.DATA_2  physical map data frame
physical_map <- read.csv(map_file, row.names = NULL, sep = ",", strip.white = TRUE)

### find markers in apr. edit physical map to include missing and remove extra markers from physical_map
lt <- list()
for(n in names(apr)){
    apr_markers <- as.data.frame(dimnames(apr[[n]])[[3]])
    colnames(apr_markers) <- c("marker")
    ##apr contains some markers not found on pmap. join on apr_markers
    physical_map_chrm <- left_join(apr_markers, physical_map, by = c("marker"))
    ##left join leaves NAs to fill
    physical_map_chrm$chr[is.na(physical_map_chrm$chr)] <- n ##fill chrm with n
    physical_map_chrm$pos[is.na(physical_map_chrm$pos)] <- 0.0001 ##fill chrm with n
    ##store new dataframe
    lt[[n]] <- physical_map_chrm
    #print(paste(n, length(tmp)))
    }
physical_map_edit <- do.call(rbind, lt)
print('physical_map')
print(dim(physical_map))
print('physical_map_edit')
print(dim(physical_map_edit))

# Convert map data frame to list
pmap <- physical_map_edit %>% 
        mutate(chr = factor(chr, levels = c(1:19, 'X'))) %>% 
        arrange(chr, pos) %>% 
        # Make sure the chromosoms are factored properly, 
        data.frame() %>% map_df_to_list(marker_column = 'marker' , pos_column = 'pos', chr_column= 'chr')

#convert pmap value from char to numeric. this is imp to use to keep markers as row names 
for(n in names(pmap)){
  pmap[[n]] <- as.list(as.data.frame(lapply(pmap[[n]], function(x) as.numeric(x))))  
  
  #make sure position is in mega_bases - make it compatibile with snp_scan --- cc_variants.sqlite data
  #pmap[[n]] <- lapply(pmap[[n]], function(x) x /1000000)
  #flatten pmap[[n]] from list to vector of numerals
  pmap[[n]] <- unlist(pmap[[n]])
}
print('pmap')
print(class(pmap[[n]] ))
print(pmap[[n]][1:3])

####################################################################
## LOAD.DATA_3 Load Phenotype Data
#pheno_data_file <- user argument 'advia_results_rankZ_20180823_kwEDIT.csv' 
#focal_phenotype <- user argument

print('load pheno_data')
pheno_data <- read.csv(pheno_data_file, sep = ",", strip.white = TRUE)
#print(head(pheno_data))
       
### create pheno_data data frame
pheno_data_df <- pheno_data %>%
  select(focal_phenotype, MouseID) 
pheno_data_df <- pheno_data_df %>%
    na.omit()   
                                            
pheno_data_df <- pheno_data_df %>% mutate(Diet=gsub("-....$","",MouseID)) 
pheno_data_df <- pheno_data_df %>% mutate(Diet=gsub("DO-","",Diet))                                       
rownames(pheno_data_df) <- pheno_data_df$MouseID    
                                            
### subset apr to animals with phenotype 
list_sample_filter <- pheno_data_df$MouseID
apr_subset <- apr[list_sample_filter]
                                     
### filter samples from pheno data -
apr_mouse_ID <- as.data.frame(dimnames(apr_subset$`1`)[[1]])
colnames(apr_mouse_ID) <- c('MouseID')
pheno_data_df <- merge(pheno_data_df, apr_mouse_ID, by = 'MouseID')
print("focal_phenotype df ")
print(dim(pheno_data_df))                                            
print(head(pheno_data_df))
                                            
####################################################################
## LOAD.DATA_4  Load Meta Data for these Mice
meta <- read.csv(meta_data_file, sep = ",", strip.white = TRUE)
print('load diet and generation wave data')
generation_tbl <- meta %>% 
  select(sample_id, gen) %>% ##edit from 191122
  rename(Generation = gen, MouseID = sample_id)
generation_tbl <- generation_tbl[!duplicated(generation_tbl),]

print('merge phenotype with generation table')
print(head(generation_tbl))
print(paste("length generation_tbl$MouseID", length(generation_tbl$MouseID)))
print(paste("length pheno_data_df$MouseID", length(pheno_data_df$MouseID))) 
                                            
pheno_data_df <- generation_tbl %>% 
  inner_join(pheno_data_df, by = 'MouseID')
rownames(pheno_data_df) <- pheno_data_df$MouseID 
print(paste("length_post_join pheno_data_df$MouseID", length(pheno_data_df$MouseID)))
print(head(pheno_data_df))

pheno_data_df$Diet <- factor(droplevels(factor(pheno_data_df$Diet)))  #use factor here to reset lelves in df
pheno_data_df$Generation <- factor(droplevels(factor(pheno_data_df$Generation)))  #use factor here to reset lelves in df
                                            
#create pheno_data for qtl
pheno_data_qtl <- pheno_data_df %>%
  select(focal_phenotype) 
rownames(pheno_data_qtl) <- pheno_data_df$MouseID  
                                            
                                            
#### 4B - Create covariate matrix
cm_DietGenwave <- model.matrix(~ Diet + Generation, data = pheno_data_df)
rownames(cm_DietGenwave) <- pheno_data_df$MouseID
cm_DietGenwave_stan <- as.matrix(cm_DietGenwave)
cm_DietGenwave <- as.matrix(cm_DietGenwave[,-1])
print(paste("dim(cm_DietGenwave)", dim(cm_DietGenwave)))

#put all cov into list , this will be accessible to focal_additive_cov
cov_list <- list(cm_DietGenwave)
names(cov_list) <- c('cm_DietGenwave')
print(paste("cov_list", names(cov_list)))
                      
####################################################################
## Analysis #1  Estiamte Kinship - run with apr not apr subset
kinship_loco <- calc_kinship(apr_subset, type = 'loco', omit_x = FALSE, cores = 1)

                                    
####################################################################
## Analysis #2  Estiamte h2

for(c in names(pmap)){  #for all chromosomes
    herit_kinbase <- est_herit(pheno_data_qtl, kinship_loco[[c]], addcovar = cov_list[[focal_additive_cov]], 
                         weights = NULL, reml = FALSE, cores = 1)
    print(paste('h2:', herit_kinbase, 'loco (leave one chrm out)',c))
} 

                                                    
                                            
####################################################################
## Analysis #2  Run additive qtl analysis                                  
genome_scan_output_file <- paste(output_base,'.genome_wide.qtl.Rdata',sep='')
print(paste("run additive qtl model. covar", focal_additive_cov))

snp_scan_data <- scan1(genoprobs = apr_subset, pheno = pheno_data_qtl, kinship = kinship_loco,
                       addcovar = cov_list[[focal_additive_cov]], cores = 1)

# Save qtl scans
save(list = c('snp_scan_data'), file = genome_scan_output_file)
                                            
#plot png  
png(paste(genome_scan_output_file,'.png',sep=''))
plot_scan1(snp_scan_data, map = pmap, type='p', cex=0.3, col=c("#775A99"), 
         main = paste0(genome_scan_output_file ))
dev.off()
        
##print to csv
qtl_csv <- paste(output_base,'.genome_wide.qtl.csv',sep='')
df_ieffects = as.data.frame(snp_scan_data)
df_ieffects$marker_name <- rownames(df_ieffects)
write.table(df_ieffects, file = qtl_csv, quote = FALSE, row.names = FALSE, sep = ',', append = FALSE)


####################################################################
####################################################################
## DATA.ANALYSIS_2B  run perumutation analysis, find signif   

alpha = 0.05
perm_thres = 1-alpha   
permute_scan_file <- paste(output_base, '.permute.',n_reps,'.Rdata',sep='')

print('run additive permutation snp scan: ')
permute_scan <- scan1perm(genoprobs = apr_subset, pheno = pheno_data_qtl, kinship = kinship_loco, 
                          addcovar = cov_list[[focal_additive_cov]], n_perm=n_reps, cores = 0)
lod_permute <- round(quantile(permute_scan,perm_thres,na.rm = TRUE), 3)
lod_permute_se <- round( (sd(permute_scan) / length(permute_scan)), 3)
print(paste("additive model", "permute reps",n_reps, "lod:", lod_permute, "se: ", lod_permute_se))

# Save permute_scan data
save(list = c('permute_scan'), file = permute_scan_file)

####################################################################
####################################################################
### DATA.ANALYSIS_3 print sig loci
sig_loci_additive <- list()
additive_threshold = lod_permute
## print significant loci
sig_loci_additive <-find_peaks(snp_scan_data, map = pmap, threshold = lod_permute, peakdrop = 1.5,drop=1.5)
print(sig_loci_additive)
## print mariginal loci
sig_loci_additive_90 <-find_peaks(snp_scan_data, map = pmap, threshold = 0.9*lod_permute, peakdrop = 1.5,drop=1.5)
print(sig_loci_additive_90)

####################################################################  
#### DATA.ANALYSIS_4 - estimate genotype effect for all 8 parental alleles    

allele_effects_file <- paste(output_base,'.genome_wide.allele_effects.Rdata',sep='')
allele_effects <- list()
print('calc linear model. parental allele_effects, additive model.')
for(c in names(pmap)){  #for all chromosomes
    print(c)

    ## rqtl2 function for BLUP - scan1blup; for linear model - scan1coef
    allele_effects[[c]] <- scan1coef(genoprobs = apr_subset[,c],    ##do not do BLUP, need to comp 
                               pheno = pheno_data_qtl,
                               kinship = kinship_loco[[c]],
                               addcovar = cov_list[[focal_additive_cov]],
                               se=TRUE,
                               cores = 1)
    }
##save allele effects
save(list = c('allele_effects'), file = allele_effects_file)
                                                                       
                                            
####################################################################
### DATA.ANALYSIS_5 -- fine mapping using imputed snps

snpscan_effects_file <- paste(output_base,'.sig_loci.snpscan_effects.Rdata',sep='')
query_variants <- create_variant_query_func(dbfile = query_variants_file)
query_genes <- create_gene_query_func(dbfile = query_genes_file)
Mb_margin <- 5  #set fine mapping window - determines the number of imputed SNPs to test. 

snp_scan_data <- list()

print('calc snp_scan additive model')
i <-  1
for(chrm in sig_loci_additive_90$chr %>% unique() %>% as.character()){

    tmp <- sig_loci_additive_90 %>% filter(chr == chrm)
    print(chrm)
    print(tmp)

    for(peak_pos in tmp$pos %>% unique() ) {

        #peak_pos = as.numeric(sig_loci_additive_90$pos[i])
        print(paste("chrm", chrm, "pos", peak_pos))
        snp_scan_data[[chrm]] <- scan1snps(genoprobs = apr_subset,
                              map = pmap,
                              pheno = pheno_data_qtl,
                              kinship = kinship_loco[[chrm]],
                              addcovar = cov_list[[focal_additive_cov]],
                              query_func = query_variants,
                              chr = chrm,
                              start = peak_pos - Mb_margin,
                              end = peak_pos + Mb_margin,
                              keep_all_snps = TRUE)
        i <-  i + 1
        ##print to tsv
        df_fine_lod <- as.data.frame(snp_scan_data[[chrm]]$lod)
        df_fine_lod$snp_id <- rownames(df_fine_lod)
        print(head(df_fine_lod))
        df_fine_snp <- as.data.frame(snp_scan_data[[chrm]]$snpinfo)
        print(head(df_fine_snp))

        df_fine_merge <- merge(df_fine_lod, df_fine_snp, by=c('snp_id'))
        print('head(df_fine_merge)')
        print(head(df_fine_merge)) 
        fine_map_file <- paste(output_base, '_chrm',chrm,'pos',peak_pos,'.fine_map.lod.tsv',sep='')
        write.table(df_fine_merge, file = fine_map_file, quote = FALSE, row.names = FALSE, sep = "\t", append = FALSE)
    }
}

##save snp_scan_data
save(list = c('snp_scan_data'), file = snpscan_effects_file)
       