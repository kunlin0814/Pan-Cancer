library(data.table)
library(tidyverse)
library(readxl)
source("/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/Burair_filtering3/Mutect1"
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/Burair_filtering3/Mutect1"
seperator <- "/"

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
retro_gene_list <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Retro_gene_finding/RetroGeneList/new_retro_gene_list_CanFam3.1.99gtf.txt",
                         header = F)

whole_wes_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt") 
exclude <- unique(unlist(whole_wes_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))
=======
=======
>>>>>>> parent of f1292fb... Update find_significat_gene.R
variant <- fread(paste(base_dir,"significant","variant_samplewise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))
gene <- fread(paste(base_dir,"significant","gene_samplewise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))
# why dup?
>>>>>>> parent of f1292fb... Update find_significat_gene.R

## check duplicated
# a <- gene[sample_names=="CCB040105"& ensembl_id=="ENSCAFG00000001781"]

# target <- mutect_after_vaf[sample_names=="CCB040105" & ensembl_id=="ENSCAFG00000001781"]
# target_combine <- target[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
# others <- mutect_after_vaf[sample_names=="CCB040105" & ensembl_id!="ENSCAFG00000001781", .(tRef,tAlt)]
# other_combine <- others[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
# testor=rbindlist(list(target_combine,other_combine))
# p_value <- fisher.test(testor)$p.value
# info <- info[,p_value:=p_value]
# info_sum <- rbindlist(list(info_sum,info))
##


# 
# Breed_info <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan-Cancer-Manuscript/Figure1/WES_WGS_merge.xlsx",
#                          sheet ="WES_WGS")
# 
# Breed_info <- setDT(Breed_info)
<<<<<<< HEAD
=======
sample_variant <- fread(paste(base_dir,"significant","variant_samplewise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))
sample_gene <- fread(paste(base_dir,"significant","gene_samplewise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))

  


Breed_info <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan-Cancer-Manuscript/Figure1/WES_WGS_merge.xlsx",
                         sheet ="WES_WGS")

Breed_info <- setDT(Breed_info)
>>>>>>> parent of 37bf719... Update find_significat_gene.R
mutect_after_vaf <- fread(paste(base_dir,"total_final_Filtering3_VAF_Mutect_withBreeds_orientBiasShaying.gz",sep =seperator))
#mutect_after_vaf <- mutect_after_vaf[,chrom_loc:=paste(chrom,pos,sep = "_")]
total_sample <- unique(mutect_after_vaf$sample_names)

=======
mutect_after_vaf <- fread(paste(base_dir,"total_final_Filtering3_VAF_Mutect_withBreeds_orientBiasShaying.gz",sep =seperator))
#mutect_after_vaf <- mutect_after_vaf[,chrom_loc:=paste(chrom,pos,sep = "_")]
total_sample <- unique(mutect_after_vaf$sample_names)

>>>>>>> parent of f1292fb... Update find_significat_gene.R
# breed <- sapply(total_sample,FUN = match_table, column="Breeds",table=Breed_info)
# mutect_after_vaf$Breeds <- breed

### samplewise variants ##
total_info_sum <- NULL
for (sample in total_sample) {
  info_sum <- NULL
  variant_loc <- mutect_after_vaf[sample_names==sample,.(chrom_loc)]
  for (i in unique(variant_loc[["chrom_loc"]]) ){
    info <- mutect_after_vaf[sample_names==sample & chrom_loc==i,]
    target <- mutect_after_vaf[sample_names==sample & chrom_loc==i, .(tRef,tAlt)]
    target_combine <- target[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    others <- mutect_after_vaf[sample_names==sample & chrom_loc!=i, .(tRef,tAlt)]
    other_combine <- others[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    testor=rbindlist(list(target_combine,other_combine))
    p_value <- fisher.test(testor)$p.value
    info <- info[,p_value:=p_value]
    info_sum <- rbindlist(list(info_sum,info))
  }
  info_sum <- info_sum[order(p_value)]
  info_sum$BH_pvalue = p.adjust(info_sum$p_value, method = "BH")
  total_info_sum <- rbindlist(list(total_info_sum,info_sum))
}
fwrite(total_info_sum,
       file = paste(base_dir,"variant_samplewise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator)
       ,col.names = T,row.names = F,
       quote = F,
       compress = "gzip",
       sep ="\t")
### samplewise ensembl_id ##

total_info_sum <- NULL
for (index in 1:length(total_sample)) {
  print(paste("processing the",index,"sample, with total samples",length(total_sample),sep = " " ))
  sample <- total_sample[index]
  info_sum <- NULL
  variant_loc <- mutect_after_vaf[sample_names==sample,.(ensembl_id)]
<<<<<<< HEAD

  each_sample_ensmbl_id <- sort(unique(variant_loc[["ensembl_id"]]))
  for (i in each_sample_ensmbl_id ){
    gene_name <- mutect_after_vaf[sample_names==sample & ensembl_id==i, .(gene_name)]$gene_name[1]

  for (i in sort(unique(variant_loc[["ensembl_id"]])) ){
    info <- mutect_after_vaf[sample_names==sample & ensembl_id==i,]

=======
  for (i in sort(unique(variant_loc[["ensembl_id"]])) ){
    info <- mutect_after_vaf[sample_names==sample & ensembl_id==i,]
>>>>>>> parent of 37bf719... Update find_significat_gene.R
    target <- mutect_after_vaf[sample_names==sample & ensembl_id==i, .(tRef,tAlt)]
    target_combine <- target[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    others <- mutect_after_vaf[sample_names==sample & ensembl_id!=i, .(tRef,tAlt)]
    other_combine <- others[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    testor=rbindlist(list(target_combine,other_combine))
    p_value <- fisher.test(testor)$p.value
    info <- info[,p_value:=p_value]
    info_sum <- rbindlist(list(info_sum,info))
  }
  info_sum <- info_sum[order(p_value)]
  info_sum$BH_pvalue = p.adjust(info_sum$p_value, method = "BH")
  total_info_sum <- rbindlist(list(total_info_sum,info_sum))
}

fwrite(total_info_sum,
       file = paste(base_dir,"gene_samplewise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator)
       ,col.names = T,row.names = F,
       quote = F,
       compress = "gzip",
       sep ="\t")

### Tumor wise
tumor_type <- unique(mutect_after_vaf$tumor_type)
### tumor_wise variants ##

total_info_sum <- NULL
for (index in 1:length(tumor_type)) {
  print(paste("processing the",index,"tumor, with total tumor",length(tumor_type),sep = " " ))
  tumor <- tumor_type[index]
  info_sum <- NULL
  variant_loc <- mutect_after_vaf[tumor_type==tumor,.(chrom_loc)]
  for (i in unique(variant_loc[["chrom_loc"]]) ){
    gene_name <- mutect_after_vaf[tumor_type==tumor & chrom_loc==i, .(gene_name)]$gene_name[1]
    ensembl_id <- mutect_after_vaf[tumor_type==tumor & chrom_loc==i, .(ensembl_id)]$ensembl_id[1]
    target <- mutect_after_vaf[tumor_type==tumor & chrom_loc==i, .(tRef,tAlt)]
    target_combine <- target[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    others <- mutect_after_vaf[tumor_type==tumor & chrom_loc!=i, .(tRef,tAlt)]
    other_combine <- others[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    testor=rbindlist(list(target_combine,other_combine))
    p_value <- fisher.test(testor)$p.value
    info <- data.table(tumor_type = tumor, 
                       chrom_loc = i, 
                       gene_name= gene_name, 
                       ensembl_id = ensembl_id,
                       p_value = p_value)
    
    info_sum <- rbindlist(list(info_sum,info))
  }
  info_sum <- info_sum[order(p_value)]
  info_sum$BH_pvalue = p.adjust(info_sum$p_value, method = "BH")
  total_info_sum <- rbindlist(list(total_info_sum,info_sum))
}


fwrite(total_info_sum,
       file = paste(base_dir,"variant_tumorwise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator)
       ,col.names = T,row.names = F,
       quote = F,
       compress = "gzip",
       sep ="\t")


### tumor_wise ensembl_id ##

### Tumor wise
tumor_type <- unique(mutect_after_vaf$tumor_type)
### tumor_wise ensmbl_id ##

total_info_sum <- NULL
for (index in 1:length(tumor_type)) {
  print(paste("processing the",index,"tumor, with total tumor",length(tumor_type),sep = " " ))
  tumor <- tumor_type[index]
  info_sum <- NULL
  variant_loc <- mutect_after_vaf[tumor_type==tumor,.(ensembl_id)]
<<<<<<< HEAD
  each_sample_ensmbl_id <- sort(unique(variant_loc[["ensembl_id"]]))
  for (i in each_sample_ensmbl_id) {

  for (i in sort(unique(variant_loc[["ensembl_id"]])) ){
    ensembl_id <- i

    gene_name <- mutect_after_vaf[tumor_type==tumor & ensembl_id==i, .(gene_name)]$gene_name[1]
=======
  for (i in sort(unique(variant_loc[["ensembl_id"]])) ){
    info <- mutect_after_vaf[tumor_type==tumor & ensembl_id==i,]
>>>>>>> parent of 37bf719... Update find_significat_gene.R
    target <- mutect_after_vaf[tumor_type==tumor & ensembl_id==i, .(tRef,tAlt)]
    target_combine <- target[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    others <- mutect_after_vaf[tumor_type==tumor & ensembl_id!=i, .(tRef,tAlt)]
    other_combine <- others[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    testor=rbindlist(list(target_combine,other_combine))
    p_value <- fisher.test(testor)$p.value
    info <- data.table(tumor_type = tumor, 
                       gene_name= gene_name, 
                       ensembl_id = ensembl_id,
                       p_value = p_value)
    
    info_sum <- rbindlist(list(info_sum,info))
  }
  info_sum <- info_sum[order(p_value)]
  info_sum$BH_pvalue = p.adjust(info_sum$p_value, method = "BH")
  total_info_sum <- rbindlist(list(total_info_sum,info_sum))
}




fwrite(total_info_sum,
       file = paste(base_dir,"gene_tumorwise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator)
      ,col.names = T,row.names = F,
       quote = F,
      compress = "gzip",
      sep ="\t")

