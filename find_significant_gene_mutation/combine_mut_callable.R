library(data.table)
library(tidyverse)
library(readxl)
# file <- fread("C:/Users/abc73/Desktop/TestTMB/CMT-2_3_java_test.txt")
remove <- c("Ensemble_id","Whole Exome")
seperator <- "/"
base_dir <- "C:/Users/abc73/Desktop/TestTMB/Mut_output"

pan_tumor_callable <- NULL
pan_tumor_genomewide <- NULL
pan_tumor_mutation_list <- NULL

tumor_type <- list.files(base_dir)

for (each_tumor in tumor_type){
  print(each_tumor)
each_tumor_file <- list.files(paste(base_dir,each_tumor,sep = seperator))

callable_mutation_inlist <- fread(paste(base_dir,each_tumor,each_tumor_file[1],sep = "/"))
callable_mutation_inlist<- callable_mutation_inlist[,(remove):=NULL]
mutation_genomewide <- fread(paste(base_dir,each_tumor,each_tumor_file[2],sep = "/"))
mutation_inlist <- fread(paste(base_dir,each_tumor,each_tumor_file[3],sep = "/"))
  
## mutation_inlist need colnames
col_names <- colnames(callable_mutation_inlist) 
colnames(mutation_inlist) <- col_names
pan_tumor_callable <- rbindlist(list(pan_tumor_callable,callable_mutation_inlist))
pan_tumor_genomewide <- rbindlist(list(pan_tumor_genomewide,mutation_genomewide))
pan_tumor_mutation_list <- rbindlist(list(pan_tumor_mutation_list,mutation_inlist))
}
## change the format, so they all use sample_names as id
colnames(pan_tumor_callable)[1] <- "sample_names"
colnames(pan_tumor_mutation_list)[1] <- "sample_names"
colnames(pan_tumor_genomewide)[1] <- "sample_names"

### append callable, genome wilde and muation info

source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/Burair_filtering3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1/"
seperator <- "/"

#retro_gene_list <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Retro_gene_finding/RetroGeneList/new_retro_gene_list_CanFam3.1.99gtf.txt",
#                         header = F)
whole_wes_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt") 
#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt")       

exclude <- unique(unlist(whole_wes_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

mutect_after_vaf <- fread(paste(base_dir,"01_29","Burair_WithBreeds_QCpass_filtering3_mutect_after_vaf.txt",sep =seperator))
mutect_after_vaf <- mutect_after_vaf[tumor_type!="UCL",]

total_ensembl_callable <- numeric(nrow(mutect_after_vaf))
total_ensembl_mut_number <- numeric(nrow(mutect_after_vaf))
total_sample_genome_wide_mut_number <- numeric(nrow(mutect_after_vaf))
total_sample_genome_wide_mut_callable <- numeric(nrow(mutect_after_vaf))

for (i in 1:nrow(mutect_after_vaf)){
  each_sample <- mutect_after_vaf[i,.(sample_names)]$sample_names
  each_ensembl <- mutect_after_vaf[i,.(ensembl_id)]$ensembl_id
  
  total_ensembl_callable[i]  <- match_table(each_sample,each_ensembl,pan_tumor_callable)
  total_ensembl_mut_number[i] <- match_table(each_sample,each_ensembl,pan_tumor_mutation_list)
  total_sample_genome_wide_mut_number[i] <- match_table(each_sample,"non_retro_PASS",pan_tumor_genomewide)
  total_sample_genome_wide_mut_callable[i] <- match_table(each_sample,"nonr_retro_callable",pan_tumor_genomewide)
  
}
mutect_after_vaf$ensembl_mut_numer <- total_ensembl_mut_number
mutect_after_vaf$ensembl_callable <- total_ensembl_callable
mutect_after_vaf$sample_genome_wide_mut_number <- total_sample_genome_wide_mut_number
mutect_after_vaf$sample_genome_wide_mut_callable <- total_sample_genome_wide_mut_callable

# write.table(mutect_after_vaf, file = "C:/Users/abc73/Desktop/mutect_noucl_vaf_withBreeds_callable.txt",
#              sep ="\t", col.names = T, row.names =F, quote = F )
# 
# fwrite(pan_tumor_callable, file = paste(base_dir,"pan_tumor_callable_table.txt",sep = seperator),
#        col.names = T, row.names = F,quote = F, sep = "\t")
# 
# 
# fwrite(pan_tumor_genomewide, file = paste(base_dir,"pan_tumor_genomewide.txt",sep = seperator),
#        col.names = T, row.names = F,quote = F, sep = "\t")
# 
# fwrite(pan_tumor_mutation_list, file = paste(base_dir,"pan_tumor_mutation_list.txt",sep = seperator),
#        col.names = T, row.names = F,quote = F, sep = "\t")





