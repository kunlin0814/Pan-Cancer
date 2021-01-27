library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
  #"/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/Burair_filtering3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/Burair_filtering3/Mutect1"
seperator <- "/"

#retro_gene_list <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Retro_gene_finding/RetroGeneList/new_retro_gene_list_CanFam3.1.99gtf.txt",
#                         header = F)

whole_wes_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt") 
             #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt")       

exclude <- unique(unlist(whole_wes_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))


## check duplicated
# a <- gene[sample_names=="CCB040105"& ensembl_id=="ENSCAFG00000001781"]

# 
# Breed_info <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan-Cancer-Manuscript/Figure1/WES_WGS_merge.xlsx",
#                          sheet ="WES_WGS")
# 
# Breed_info <- setDT(Breed_info)
mutect_after_vaf <- fread(paste(base_dir,"total_final_Filtering3_VAF_Mutect_withBreeds_orientBiasShaying.gz",sep =seperator))

mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude & gene_name!="-", ]
total_sample <- unique(mutect_after_vaf$sample_names)


# breed <- sapply(total_sample,FUN = match_table, column="Breeds",table=Breed_info)
# mutect_after_vaf$Breeds <- breed

### samplewise variants ##
total_info_sum <- NULL
for (index in 1:length(total_sample)) {
  print(paste("processing the",index,"sample, with total samples",length(total_sample),sep = " " ))
  sample <- total_sample[index]
  info_sum <- NULL
  variant_loc <- mutect_after_vaf[sample_names==sample,.(chrom_loc)]
  for (i in unique(variant_loc[["chrom_loc"]]) ){
    info <- mutect_after_vaf[sample_names==sample & chrom_loc==i,]
    target <- mutect_after_vaf[sample_names==sample & chrom_loc==i, .(tRef,tAlt)]
    target_combine <- target[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    others <- mutect_after_vaf[sample_names==sample & chrom_loc!=i, .(tRef,tAlt)]
    other_combine <- others[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    testor=rbindlist(list(target_combine,other_combine))
    p_value <- fisher.test(testor,alternative = "less")$p.value
    info <- info[,p_value:=p_value]
    info_sum <- rbindlist(list(info_sum,info))
  }
  info_sum <- info_sum[order(p_value)]
  info_sum$BH_pvalue = p.adjust(info_sum$p_value, method = "BH")
  total_info_sum <- rbindlist(list(total_info_sum,info_sum))
}

# fwrite(total_info_sum,
#        file = paste(base_dir,"variant_samplewise_p_value_Filtering3_VAF_Mutect_orientBias3_01_26.gz",sep = seperator)
#        ,col.names = T,
#        row.names = F,
#        quote = F,
#        eol = "\n",
#        compress = "gzip",
#        sep ="\t")
### samplewise ensembl_id ##


gene_total_info_sum <- NULL
for (index in 1:length(total_sample)) {
  print(paste("processing the",index,"sample, with total samples",length(total_sample),sep = " " ))
  sample <- total_sample[index]
  info_sum <- NULL
  tumor <-  mutect_after_vaf[sample_names==sample,.(tumor_type)]$tumor_type[1]
  ensembl_name <- mutect_after_vaf[sample_names==sample,.(ensembl_id)]
  each_sample_ensmbl_id <- sort(unique(ensembl_name[["ensembl_id"]]))
  for (i in each_sample_ensmbl_id ){
    gene_name <- mutect_after_vaf[sample_names==sample & ensembl_id==i, .(gene_name)]$gene_name[1]
    target <- mutect_after_vaf[sample_names==sample & ensembl_id==i, .(tRef,tAlt)]
    target_combine <- target[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    others <- mutect_after_vaf[sample_names==sample & ensembl_id!=i, .(tRef,tAlt)]
    other_combine <- others[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    testor=rbindlist(list(target_combine,other_combine))
    p_value <- fisher.test(testor,alternative = "less")$p.value
    info <- data.table(sample_names = sample, 
                       gene_name= gene_name, 
                       ensembl_id = i,
                       tumor_type = tumor,
                       p_value = p_value)
    info_sum <- rbindlist(list(info_sum,info))
  }
  info_sum <- info_sum[order(p_value)]
  info_sum$BH_pvalue = p.adjust(info_sum$p_value, method = "BH")
  gene_total_info_sum <- rbindlist(list(gene_total_info_sum,info_sum))
}

# fwrite(gene_total_info_sum,
#        file = paste(base_dir,"gene_samplewise_p_value_Filtering3_VAF_Mutect_orientBias3_01_26.gz",sep = seperator)
#        ,col.names = T,row.names = F,
#        quote = F,
#        eol = "\n",
#        compress = "gzip",
#        sep ="\t")

## TumorWise
tumor_type <- unique(mutect_after_vaf$tumor_type)
### Tumorwise variants

sig_variants <- total_info_sum[BH_pvalue < 0.05,]
total_variant_summary <- NULL
for (index in 1:length(tumor_type)) {
  #print(paste("processing the", index, "tumor, with total number tumor", length(tumor_type),sep = " "))
  each_tumor <- tumor_type[index]
  each_tumor_info_sum <- NULL
  each_tumor_type <-
    sig_variants[tumor_type == each_tumor, .(chrom_loc, sample_names, gene_name, ensembl_id)]
  total_variants_each_tumor <- unique(each_tumor_type[["chrom_loc"]])
  for (i in 1:length(total_variants_each_tumor)) {
    #print(paste("processing the",i, "variants, with total variants", length(total_variants_each_tumor), sep = " "))
    each_variant <- total_variants_each_tumor[i]
    total_others <- 0
    total_others_without <- 0
    numbersamples_withvariants <-
      length(unique(each_tumor_type[chrom_loc == each_variant, .(sample_names)][['sample_names']]))
    numersamples_withoutvariants <-
      length(unique(each_tumor_type$sample_names)) - numbersamples_withvariants
    gene_name <-
      each_tumor_type[chrom_loc == each_variant, .(gene_name)][['gene_name']][1]
    ensembl_id <-
      each_tumor_type[chrom_loc == each_variant, .(ensembl_id)][['ensembl_id']][1]
    
    remaining_variant <- total_variants_each_tumor[!total_variants_each_tumor %in% each_variant]
    
    for (other_variant in  remaining_variant) {
    
        nubmersamples_others_variants <- length(unique(each_tumor_type[chrom_loc == other_variant, .(sample_names)][['sample_names']]))
        numersamples_others_withoutvariants <- length(unique(each_tumor_type$sample_names)) - nubmersamples_others_variants
        total_others <- total_others + nubmersamples_others_variants
        total_others_without <- total_others_without + numersamples_others_withoutvariants
      
      
    }
    testor <-
      rbind(c(numbersamples_withvariants,numersamples_withoutvariants),
        c(total_others, total_others_without))
    
    each_variant_p_value <-
      fisher.test(testor, alternative = "greater")$p.value
    # if (each_variant_p_value <0.05){
    # print(testor)
    # print(each_variant_p_value)}
    info <- data.table(
      tumor_type = each_tumor,
      chrom_loc = each_variant,
      gene_name = gene_name,
      ensembl_id = ensembl_id,
      p_value = each_variant_p_value
    )
    
    each_tumor_info_sum <- rbindlist(list(each_tumor_info_sum, info))
  }
  each_tumor_info_sum <- each_tumor_info_sum[order(p_value)]
  each_tumor_info_sum$BH_pvalue = p.adjust(each_tumor_info_sum$p_value, method = "BH")
  total_variant_summary <- rbindlist(list(total_variant_summary, each_tumor_info_sum))
}

fwrite(total_variant_summary,
       file = paste(base_dir,"variant_tumorwise_p_value_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator)
       ,col.names = T,row.names = F,
       eol = "\n",
       quote = F,
       compress = "gzip",
       sep ="\t")

a <- total_variant_summary[tumor_type=="MT"]

### Tumorwise genes
sig_gene <- gene_total_info_sum[BH_pvalue<0.05,]

total_ensembl_summary <- NULL
for (index in 1:length(tumor_type)) {
  print(paste("processing the", index, "tumor, with total number tumor", length(tumor_type),sep = " "))
  each_tumor <- tumor_type[index]
  each_tumor_info_sum <- NULL
  each_tumor_type <- sig_gene[tumor_type == each_tumor, .(sample_names, gene_name, ensembl_id)]
  total_ensembl_each_tumor <- unique(each_tumor_type[["ensembl_id"]])
  for (i in 1:length(total_ensembl_each_tumor)) {
    print(paste("processing the",i, "ensembl, with total ensembl", length(total_ensembl_each_tumor), sep = " "))
    each_ensmbl <- total_ensembl_each_tumor[i]
    total_others <- 0
    total_others_without <- 0
    numbersamples_withvariants <-
      length(unique(each_tumor_type[ensembl_id == each_ensmbl, .(sample_names)][['sample_names']]))
    numersamples_withoutvariants <-
      length(unique(each_tumor_type$sample_names)) - numbersamples_withvariants
    gene_name <- each_tumor_type[ensembl_id == each_ensmbl, .(gene_name)][['gene_name']][1]
    
    for (other_ensembl in  total_ensembl_each_tumor) {
      if (other_ensembl != each_ensmbl) {
        nubmersamples_others_variants <- length(unique(each_tumor_type[ensembl_id == other_ensembl, .(sample_names)][['sample_names']]))
        numersamples_others_withoutvariants <- length(unique(each_tumor_type$sample_names)) - nubmersamples_others_variants
        
        total_others <- total_others + nubmersamples_others_variants
        total_others_without <- total_others_without + numersamples_others_withoutvariants
      }
      
    }
    testor <-
      rbind(c(numbersamples_withvariants,numersamples_withoutvariants),
            c(total_others, total_others_without))
    
    each_ensembl_p_value <-
      fisher.test(testor, alternative = "greater")$p.value
    
    info <- data.table(
      tumor_type = each_tumor,
      gene_name = gene_name,
      ensembl_id = each_ensmbl,
      p_value = each_ensembl_p_value
    )
    
    each_tumor_info_sum <- rbindlist(list(each_tumor_info_sum, info))
  }
  each_tumor_info_sum <- each_tumor_info_sum[order(p_value)]
  each_tumor_info_sum$BH_pvalue = p.adjust(each_tumor_info_sum$p_value, method = "BH")
  total_ensembl_summary <- rbindlist(list(total_ensembl_summary, each_tumor_info_sum))
}


# fwrite(total_ensembl_summary,
#        file = paste(base_dir,"gene_tumorwise_p_value_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator)
#       ,col.names = T,row.names = F,
#        quote = F,
#       eol = "\n",
#       compress = "gzip",
#       sep ="\t")

# 
# ### check the results
# # 
# sample_variant <- fread(paste(base_dir,"significant","clean_variant_samplewise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))
# sample_gene <- fread(paste(base_dir,"significant","clean_gene_samplewise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))
# tumor_variant <- fread(paste(base_dir,"significant","clean_variant_tumorwise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))
# tumor_gene <- fread(paste(base_dir,"significant","clean_gene_tumorwise_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))
# 
# a <- sample_variant[sample_names=="004",]
# 
# significant_sample_variant <- sample_variant[!ensembl_id %in% retro_gene_list
#                                              $V1 & BH_pvalue < 0.05 
#                                              & gene_name!="-",.(chrom_loc, BH_pvalue,vaf), 
#                                              keyby = .(sample_names,tumor_type)]
# 
# significant_sample_gene <- sample_variant[!ensembl_id %in% retro_gene_list$V1 
#                                           & BH_pvalue < 0.05
#                                           & gene_name!="-",.(gene_name, BH_pvalue), 
#                                           keyby = .(sample_names,tumor_type)]
# significant_tumor_variant <- tumor_variant[!ensembl_id %in% retro_gene_list$V1 
#                                            & BH_pvalue < 0.05
#                                            & gene_name!="-",.(chrom_loc, BH_pvalue),
#                                            keyby = .(tumor_type)]
# significant_tumor_gene <- tumor_gene[!ensembl_id %in% retro_gene_list$V1 
#                                      & BH_pvalue < 0.05
#                                      & gene_name!="-",
#                                      .(gene_name, BH_pvalue), keyby = .(tumor_type)]
# 
# 
# top_10_sample_variants <- significant_sample_variant[, head(.SD, 10), by=.(sample_names)]
# top_10_sample_genes <- significant_sample_gene[, head(.SD, 10), by=.(sample_names)]
# top_10_tumor_variants <- significant_tumor_variant[, head(.SD, 10), by=tumor_type]
# top_10_tumor_gene <- significant_tumor_gene[, head(.SD, 10), by=tumor_type]
# 
# 
# # a <- top_10_sample_variants[tumor_type=="HSA",]
# # 
# # 
