library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
seperator <- "/"

#retro_gene_list <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Retro_gene_finding/RetroGeneList/new_retro_gene_list_CanFam3.1.99gtf.txt",
#                         header = F)
whole_wes_table <- fread(#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt") 
"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt")       

exclude <- unique(unlist(whole_wes_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

amp_delete <- fread(paste(base_dir,"With_subtype_no_pesudo_samples_amp_delete_0202.txt",sep = seperator))

#colnames(amp_delete) <- c("sample_names","gene_name","mutation_type","CNA","tumor_type")
# 
# amp_delete <- amp_delete[!sample_names %in% exclude]
# amp_delete <- amp_delete[,gene_mutation:=paste(gene_name,mutation_type,sep = "_")]
# 
# ## append subtype 
# table_total_sample <- amp_delete$sample_names
# subtype <- sapply(table_total_sample,FUN = match_table, column="DiseaseAcronym2",table=whole_wes_table)
# amp_delete$subtype <- subtype
# ## get rid of pseuo genes
# amp_delete <- amp_delete[!grepl("ENSCAFG",amp_delete[,.(gene_name)]$gene_name,ignore.case = T)]
# 
# fwrite(amp_delete,file = paste(base_dir,"With_subtype_no_pesudo_samples_amp_delete_0202.txt",sep = seperator),
#        col.names = T, row.names = F, quote = F, sep = "\t",
#        eol = "\n")

################ main code ####################
tumor_type <- sort(unique(amp_delete$subtype))
## use which to find index and check how many samples ( can't use unique gene to find)
## 2/2 merge amp and delete
## 2/3 seperate merge amp and delete and treat them as different mut


total_gene_summary <- NULL
for (each_tumor in tumor_type) {
  current <- match(each_tumor, tumor_type,nomatch = 0)
  print(paste("Processing the ", current, "tumors with total numbers of tumor", length(tumor_type),sep = " "))
  each_tumor_info_sum <- NULL
  each_tumor_info <- amp_delete[subtype == each_tumor, ]
  each_tumor_mut_type <- each_tumor_info$gene_mutation
  each_tumor_uniq_mut_type <- unique(each_tumor_mut_type)
  each_tumor_total_sample <-   length(unique(each_tumor_info$sample_names))
  
  tumor_type_sum <- character(length(each_tumor_uniq_mut_type))
  gene_mut_sum <- character(length(each_tumor_uniq_mut_type))
  gene_sum <- character(length(each_tumor_uniq_mut_type))
  numbersamples_withgene_sum <-numeric(length(each_tumor_uniq_mut_type))
  numbersamples_withoutgene_sum <- numeric(length(each_tumor_uniq_mut_type))
  total_others_sum <- numeric(length(each_tumor_uniq_mut_type))
  total_others_without_sum <- numeric(length(each_tumor_uniq_mut_type))
  
  total_numbersamples_withgene <- 0
  total_numbersamples_withoutgene <- 0 
  
  # calculate total x and total y
  for (index in 1:length(each_tumor_uniq_mut_type)) {
    #current <- match(each_gene, each_tumor_gene,nomatch = 0)
    print(paste("Processing the ", index, "gene mutation, with total gene mutations", length(each_tumor_uniq_mut_type),sep = " "))
    each_gene <- each_tumor_uniq_mut_type[index]
    total_others <- 0
    total_others_without <- 0
    sample_loc <- which(each_gene == each_tumor_mut_type)
    #length(unique(each_tumor_info[gene_mutation == each_gene, ]$sample_names))
    numbersamples_withgene <- length(unique(each_tumor_info$sample_names[sample_loc]))
    #length(which(each_gene == each_tumor_gene))
    #length(unique(each_tumor_info[gene_mutation == each_gene, ]$sample_names))
    
    numbersamples_withoutgene <- each_tumor_total_sample- numbersamples_withgene
    outside_target_gene <- each_tumor_uniq_mut_type[!each_tumor_uniq_mut_type %in% each_gene]
    
    total_numbersamples_withgene <- total_numbersamples_withgene+numbersamples_withgene
    total_numbersamples_withoutgene <- total_numbersamples_withoutgene+numbersamples_withoutgene
    
    #gene_sum[index] <- each_gene
    tumor_type_sum[index] <- each_tumor
    gene_mut_sum[index] <- each_gene
    numbersamples_withgene_sum[index] <- numbersamples_withgene
    numbersamples_withoutgene_sum[index] <- numbersamples_withoutgene
    
    
  }
  total_others_sum <- total_numbersamples_withgene -numbersamples_withgene_sum
  total_others_without_sum <- total_numbersamples_withoutgene- numbersamples_withoutgene_sum
  
  
  
  
  each_tumor_info_sum <- data.table( tumor_type = tumor_type_sum,
                                     gene_mut = gene_mut_sum,
                                     numbersamples_withgene = numbersamples_withgene_sum,
                                     numbersamples_withoutgene = numbersamples_withoutgene_sum,
                                     total_others = total_others_sum,
                                     total_others_withoutgene = total_others_without_sum
  )
  
  p_value_sum <- numeric(length(each_tumor_uniq_mut_type))
  for ( i in 1:nrow(each_tumor_info_sum)) {
    
    target <- as.matrix(each_tumor_info_sum[i,.(numbersamples_withgene,numbersamples_withoutgene)])
    others <- as.matrix(each_tumor_info_sum[i,.(total_others,total_others_withoutgene)])
    testor <- rbind(c(target),c(others))
    each_gene_p_value <- fisher.test(testor, alternative = "greater")$p.value
    p_value_sum[i] <-  each_gene_p_value
    
  }
  each_tumor_info_sum$p_value <- p_value_sum 
  total_gene_summary <- rbindlist(list(total_gene_summary, each_tumor_info_sum))
}



### Analyzed the results, sep gene and mut type

amp_delete <- total_gene_summary

different_type <- as.data.table(str_split_fixed(amp_delete$gene_mut,"_",2))
colnames(different_type) <- c("gene_name","mut_type")
amp_delete <- cbind(amp_delete,different_type)

# 
# fwrite(amp_delete, file = paste(base_dir, "sep_amp_delete_before_BH_Pan_cancer_amp_delete_02_03.txt",
#                                         sep = seperator),
#        col.names = T, row.names = F, quote = F,sep = "\t",
#        eol = "\n")

#### adjust pvalue with BH methods ####


total_summary <- NULL
tumor_type <- sort(unique(amp_delete$tumor_type))

for (each_tumor in tumor_type) {
    each_tumor_sum <- NULL
    each_tumor_info <- amp_delete[tumor_type==each_tumor]
    #gene_name <- sort(unique(each_tumor_info[,.(gene_name)]$gene_name))
    #each_mut_type_tumor <- each_tumor_info[gene_name==each_mut_type,]
    each_mut_type_tumor <- each_tumor_info[order(p_value)]
    each_mut_type_tumor$BH_pvalue <- p.adjust(each_mut_type_tumor$p_value, method = "BH")
    each_tumor_sum <- rbindlist(list(each_tumor_sum,each_mut_type_tumor))
    

  total_summary <- rbindlist(list(total_summary,each_tumor_sum))
}

# fwrite(total_summary, file = paste(base_dir,"02_03", "With_BH_sep_amp_delete_merged_adjust_Pan_cancer_amp_delete_02_03.txt",
#                                         sep = seperator),
#        col.names = T, row.names = F, quote = F,sep = "\t",
#        eol = "\n")

### Analyzed the final results ###
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
base_dir <- "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
seperator <- "/"

original_amp_delete <- fread(paste(base_dir, "With_subtype_no_pesudo_samples_amp_delete_0202.txt", sep = seperator))

p_value_amp_delete <- fread(paste(base_dir,"02_03", "With_BH_sep_amp_delete_merged_adjust_Pan_cancer_amp_delete_02_03.txt", sep = seperator))                                    
merge_gene_amp_delete <- fread(paste(base_dir,"02_02", "With_BH_final_Pan_cancer_amp_delete_02_02.txt", sep = seperator))                                    

# check the result
top10_sep <- p_value_amp_delete[,head(.SD,10),keyby = tumor_type] 
top10_merge <- merge_gene_amp_delete[,head(.SD,10),keyby = tumor_type]
MT <- p_value_amp_delete[tumor_type=="MT" ]
a <- MT[grepl("CDKN",p_value_amp_delete[tumor_type=="MT",.(gene_name)][["gene_name"]], ignore.case = T)]

sig_amp_del <- p_value_amp_delete[BH_pvalue < 0.05]

# mut_type <- amp_delete[ ,.N, keyby = .(tumor_type,mut_type)]



