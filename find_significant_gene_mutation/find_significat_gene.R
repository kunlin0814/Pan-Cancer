library(data.table)
library(tidyverse)
library(readxl)
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


## check duplicated
# a <- gene[sample_names=="CCB040105"& ensembl_id=="ENSCAFG00000001781"]

# 
Breed_info <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan-Cancer-Manuscript/Figure1/WES_WGS_merge.xlsx",
                         sheet ="WES_WGS")

Breed_info <- setDT(Breed_info)
mutect_after_vaf <- fread(paste(base_dir,"01_31","mutect_noucl_vaf_withBreeds_callable.txt",sep =seperator))

# mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude & gene_name!="-", ]
# 
# table_total_sample <- mutect_after_vaf$sample_names
# 
# breed <- sapply(table_total_sample,FUN = match_table, column="Breeds",table=Breed_info)
# mutect_after_vaf$Breeds <- breed
# mutect_after_vaf <- mutect_after_vaf[,chrom_loc:= paste(chrom,pos,sep = "_"),]

# write.table(mutect_after_vaf, file = "C:/Users/abc73/Desktop/Burair_WithBreeds_QCpass_filtering3_mutect_after_vaf.txt",
#             sep ="\t", col.names = T, row.names =F, quote = F )

# a <- unique(mutect_after_vaf$ensembl_id)
# write.table(a, file = "C:/Users/abc73/Desktop/total_target_ensmbl_id",
#             sep ="\n", col.names = T, row.names =F, quote = F )


total_sample <- unique( mutect_after_vaf$sample_names)
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

fwrite(total_info_sum,
       file = paste(base_dir,"01_31","variant_samplewise_p_value_Filtering3_VAF_Mutect_orientBias3_01_31.gz",sep = seperator)
       ,col.names = T,
       row.names = F,
       quote = F,
       eol = "\n",
       compress = "gzip",
       sep ="\t")


### samplewise ensembl_id ##


gene_total_info_sum <- NULL
for (index in 1:length(total_sample)) {
  print(paste("processing the",index,"sample, with total samples",length(total_sample),sep = " " ))
  sample <- total_sample[index]
  info_sum <- NULL
  tumor <-  mutect_after_vaf[sample_names==sample,.(tumor_type)]$tumor_type[1]
  gene_name <- mutect_after_vaf[sample_names==sample,.(gene_name)]
  each_sample_gene_id <- sort(unique(gene_name[["gene_name"]]))
  for (i in each_sample_gene_id ){
    ensembl_id <- mutect_after_vaf[sample_names==sample & gene_name==i, .(ensembl_id)]$ensembl_id[1]
    target <- mutect_after_vaf[sample_names==sample & gene_name==i, .(tRef,tAlt)]
    target_combine <- target[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    others <- mutect_after_vaf[sample_names==sample & gene_name!=i, .(tRef,tAlt)]
    other_combine <- others[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    testor=rbindlist(list(target_combine,other_combine))
    p_value <- fisher.test(testor,alternative = "less")$p.value
    info <- data.table(sample_names = sample, 
                       gene_name= i, 
                       ensembl_id = ensembl_id,
                       tumor_type = tumor,
                       p_value = p_value)
    info_sum <- rbindlist(list(info_sum,info))
  }
  info_sum <- info_sum[order(p_value)]
  info_sum$BH_pvalue = p.adjust(info_sum$p_value, method = "BH")
  gene_total_info_sum <- rbindlist(list(gene_total_info_sum,info_sum))
}

fwrite(gene_total_info_sum,
       file = paste(base_dir,"01_31","pure_gene_samplewise_p_value_Filtering3_VAF_Mutect_orientBias3_01_31.gz",sep = seperator)
       ,col.names = T,row.names = F,
       quote = F,
       eol = "\n",
       compress = "gzip",
       sep ="\t")


sig_variants <- fread(file = paste(base_dir,"01_31","variant_samplewise_p_value_Filtering3_VAF_Mutect_orientBias3_01_31.gz",sep = seperator))
sig_variants <- sig_variants[BH_pvalue < 0.05,]

sig_gene <- fread(paste(base_dir,"01_31","pure_gene_samplewise_p_value_Filtering3_VAF_Mutect_orientBias3_01_31.gz",sep = seperator))
sig_gene <- sig_gene[BH_pvalue < 0.05, ]

# fwrite(sig_gene,
#        file = paste(base_dir,"gene_samplewise_pvaluelt0.05_Filtering3_VAF_Mutect_orientBias3_01_26.gz",sep = seperator)
#        ,col.names = T,row.names = F,
#        quote = F,
#        eol = "\n",
#        compress = "gzip",
#        sep ="\t")

## TumorWise
tumor_type <- unique(sig_variants$tumor_type)
### Tumorwise variants

total_variant_summary <- NULL
for (each_tumor in tumor_type) {
  current <- match(each_tumor, tumor_type,nomatch = 0)
  print(paste("Processing the ", current, "tumors with total numbers of tumor", length(tumor_type),sep = " "))
  each_tumor_info_sum <- NULL
  each_tumor_info <- sig_variants[tumor_type == each_tumor, ]
  each_tumor_variants <- each_tumor_info$chrom_loc
  each_tumor_uniq_variants <- unique(each_tumor_variants)
  tumor_type_sum <- character(length(each_tumor_uniq_variants))
  gene_mut_sum <- character(length(each_tumor_uniq_variants))
  numbersamples_withgene_sum <-numeric(length(each_tumor_uniq_variants))
  numbersamples_withoutgene_sum <- numeric(length(each_tumor_uniq_variants))
  total_others_sum <- numeric(length(each_tumor_uniq_variants))
  total_others_without_sum <- numeric(length(each_tumor_uniq_variants))
  each_tumor_total_sample <-   length(unique(each_tumor_info$sample_names))
  ensmbl_id_sum <- character(length(each_tumor_uniq_variants))
  gene_name_sum <- character(length(each_tumor_uniq_variants))
  total_numbersamples_withgene <- 0
  total_numbersamples_withoutgene <- 0 
  
  # calculate total x and total y
  for (index in 1:length(each_tumor_uniq_variants)) {
    #current <- match(each_variant, each_tumor_variants,nomatch = 0)
    print(paste("Processing the ", index, "variants, with variants", length(each_tumor_uniq_variants),sep = " "))
    each_variant <- each_tumor_uniq_variants[index]
    each_gene <- each_tumor_info[chrom_loc==each_variant]$gene_name[1]
    each_ensembl <- each_tumor_info[chrom_loc==each_variant]$ensembl_id[1]
    total_others <- 0
    total_others_without <- 0
    sample_loc <- which(each_variant == each_tumor_variants)
    
    #length(unique(each_tumor_info[gene_mutation == each_variant, ]$sample_names))
    numbersamples_withgene <- length(unique(each_tumor_info$sample_names[sample_loc]))
    #length(which(each_variant == each_tumor_variants))
    #length(unique(each_tumor_info[gene_mutation == each_variant, ]$sample_names))
    
    numbersamples_withoutgene <- each_tumor_total_sample- numbersamples_withgene
    
    total_numbersamples_withgene <- total_numbersamples_withgene+numbersamples_withgene
    total_numbersamples_withoutgene <- total_numbersamples_withoutgene+numbersamples_withoutgene
    
    
    gene_name_sum[index] <- each_gene
    ensmbl_id_sum[index] <- each_ensembl
    tumor_type_sum[index] <- each_tumor
    gene_mut_sum[index] <- each_variant
    numbersamples_withgene_sum[index] <- numbersamples_withgene
    numbersamples_withoutgene_sum[index] <- numbersamples_withoutgene
    
    
  }
  total_others_sum <- total_numbersamples_withgene -numbersamples_withgene_sum
  total_others_without_sum <- total_numbersamples_withoutgene- numbersamples_withoutgene_sum
  
  
  
  
  each_tumor_info_sum <- data.table( tumor_type = tumor_type_sum,
                                     ensembl_id =ensmbl_id_sum,
                                     gene_name = gene_name_sum,
                                     gene_mut = gene_mut_sum,
                                     numbersamples_withgene = numbersamples_withgene_sum,
                                     numbersamples_withoutgene = numbersamples_withoutgene_sum,
                                     total_others = total_others_sum,
                                     total_others_withoutgene = total_others_without_sum
  )
  
  p_value_sum <- numeric(length(each_tumor_uniq_variants))
  for ( i in 1:nrow(each_tumor_info_sum)) {
    
    target <- as.matrix(each_tumor_info_sum[i,.(numbersamples_withgene,numbersamples_withoutgene)])
    others <- as.matrix(each_tumor_info_sum[i,.(total_others,total_others_withoutgene)])
    testor <- rbind(c(target),c(others))
    each_gene_p_value <- fisher.test(testor, alternative = "greater")$p.value
    p_value_sum[i] <-  each_gene_p_value
    
  }
  each_tumor_info_sum$p_value <- p_value_sum 
  each_tumor_info_sum <- each_tumor_info_sum[order(p_value)]
  each_tumor_info_sum$BH_pvalue = p.adjust(each_tumor_info_sum$p_value, method = "BH")
  
  total_variant_summary <- rbindlist(list(total_variant_summary, each_tumor_info_sum))
}

fwrite(total_variant_summary,
       file = paste(base_dir,"01_31","variant_tumorwise_p_value_Filtering3_VAF_Mutect_orientBias3_0129.gz",sep = seperator)
       ,col.names = T,row.names = F,
       eol = "\n",
       quote = F,
       compress = "gzip",
       sep ="\t")


## TumorWise
tumor_type <- unique(sig_gene$tumor_type)
### Tumorwise genes
total_gene_summary <- NULL
for (each_tumor in tumor_type) {
  current <- match(each_tumor, tumor_type,nomatch = 0)
  print(paste("Processing the ", current, "tumors with total numbers of tumor", length(tumor_type),sep = " "))
  each_tumor_info_sum <- NULL
  each_tumor_info <- sig_gene[tumor_type == each_tumor, ]
  each_tumor_genes <- each_tumor_info$gene_name
  each_tumor_uniq_genes <- unique(each_tumor_genes)
  tumor_type_sum <- character(length(each_tumor_uniq_genes))
  numbersamples_withgene_sum <-numeric(length(each_tumor_uniq_genes))
  numbersamples_withoutgene_sum <- numeric(length(each_tumor_uniq_genes))
  total_others_sum <- numeric(length(each_tumor_uniq_genes))
  total_others_without_sum <- numeric(length(each_tumor_uniq_genes))
  each_tumor_total_sample <-   length(unique(each_tumor_info$sample_names))
  ensmbl_id_sum <- character(length(each_tumor_uniq_genes))
  gene_name_sum <- character(length(each_tumor_uniq_genes))
  total_numbersamples_withgene <- 0
  total_numbersamples_withoutgene <- 0 
  
  # calculate total x and total y
  for (index in 1:length(each_tumor_uniq_genes)) {
    #current <- match(each_gene, each_tumor_genes,nomatch = 0)
    print(paste("Processing the ", index, "genes, with total_genes", length(each_tumor_uniq_genes),sep = " "))
    each_gene <- each_tumor_uniq_genes[index]
    each_ensembl <- each_tumor_info[gene_name==each_gene]$ensembl_id[1]
    total_others <- 0
    total_others_without <- 0
    sample_loc <- which(each_gene == each_tumor_genes)
    
    #length(unique(each_tumor_info[gene_mutation == each_gene, ]$sample_names))
    numbersamples_withgene <- length(unique(each_tumor_info$sample_names[sample_loc]))
    #length(which(each_gene == each_tumor_genes))
    #length(unique(each_tumor_info[gene_mutation == each_gene, ]$sample_names))
    
    numbersamples_withoutgene <- each_tumor_total_sample- numbersamples_withgene
    
    total_numbersamples_withgene <- total_numbersamples_withgene+numbersamples_withgene
    total_numbersamples_withoutgene <- total_numbersamples_withoutgene+numbersamples_withoutgene
    
    
    gene_name_sum[index] <- each_gene
    ensmbl_id_sum[index] <- each_ensembl
    tumor_type_sum[index] <- each_tumor
    numbersamples_withgene_sum[index] <- numbersamples_withgene
    numbersamples_withoutgene_sum[index] <- numbersamples_withoutgene
    
    
  }
  total_others_sum <- total_numbersamples_withgene -numbersamples_withgene_sum
  total_others_without_sum <- total_numbersamples_withoutgene- numbersamples_withoutgene_sum
  
  
  
  
  each_tumor_info_sum <- data.table( tumor_type = tumor_type_sum,
                                     ensembl_id =ensmbl_id_sum,
                                     gene_name = gene_name_sum,
                                     numbersamples_withgene = numbersamples_withgene_sum,
                                     numbersamples_withoutgene = numbersamples_withoutgene_sum,
                                     total_others = total_others_sum,
                                     total_others_withoutgene = total_others_without_sum
  )
  
  p_value_sum <- numeric(length(each_tumor_uniq_genes))
  for ( i in 1:nrow(each_tumor_info_sum)) {
    
    target <- as.matrix(each_tumor_info_sum[i,.(numbersamples_withgene,numbersamples_withoutgene)])
    others <- as.matrix(each_tumor_info_sum[i,.(total_others,total_others_withoutgene)])
    testor <- rbind(c(target),c(others))
    each_gene_p_value <- fisher.test(testor, alternative = "greater")$p.value
    p_value_sum[i] <-  each_gene_p_value
    
  }
  each_tumor_info_sum$p_value <- p_value_sum 
  each_tumor_info_sum <- each_tumor_info_sum[order(p_value)]
  each_tumor_info_sum$BH_pvalue = p.adjust(each_tumor_info_sum$p_value, method = "BH")
  
  total_gene_summary <- rbindlist(list(total_gene_summary, each_tumor_info_sum))
}

check <- total_ensembl_summary[BH_pvalue<0.2,]



fwrite(total_ensembl_summary,
       file = paste(base_dir,"01_31","pure_gene_tumorwise_p_value_Filtering3_VAF_Mutect_orientBias3_01_27.gz",sep = seperator)
      ,col.names = T,row.names = F,
       quote = F,
      eol = "\n",
      compress = "gzip",
      sep ="\t")

# 
# ### check the results

variant_sample <- fread(paste(base_dir,"01_27","variant_samplewise_p_value_Filtering3_VAF_Mutect_orientBias3_01_27.gz",
                              sep = seperator))
variant_tumor <- fread(paste(base_dir,"01_27","variant_tumorwise_p_value_Filtering3_VAF_Mutect_orientBias3_01_27.gz",
                             sep = seperator))
gene_sample <- fread(paste(base_dir,"01_27","gene_samplewise_p_value_Filtering3_VAF_Mutect_orientBias3_01_27.gz",
                           sep = seperator))
gene_tumor <- fread(paste(base_dir,"01_27","gene_tumorwise_p_value_Filtering3_VAF_Mutect_orientBias3_01_27.gz",
                          sep = seperator))

pik3_sign <- unique(variant_sample[tumor_type=="HSA" & gene_name=="PIK3CA",.(sample_names)])
tp53_sign <- unique(variant_sample[tumor_type=="HSA" & gene_name=="TP53",.(sample_names)])

freq <- as.data.frame(table(tumor$chrom_loc))

file <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/Burair_filtering3/Mutect1/01_26/variant_tumorwise_p_value_Filtering3_VAF_Mutect_orientBias3.gz") 
tumor <- file[tumor_type=="HSA",] 

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
