library(data.table)
library(tidyverse)
library(readxl)
source(#"C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
"/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
seperator <- "/"

#retro_gene_list <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Retro_gene_finding/RetroGeneList/new_retro_gene_list_CanFam3.1.99gtf.txt",
#                         header = F)
whole_wes_table <- fread(#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt") 
"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt")       

exclude <- unique(unlist(whole_wes_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

amp_delete <- fread(paste(base_dir,"total_samples_amp_delete.txt",sep = seperator),header = F)
colnames(amp_delete) <- c("sample_names","gene_name","mutation_type","CNA","tumor_type")

amp_delete <- amp_delete[!sample_names %in% exclude,gene_mutation:=paste(gene_name,mutation_type,sep = "_")]


tumor_type <- sort(unique(amp_delete$tumor_type))
## use which to find index and check how many samples ( can't use unique gene to find)
total_gene_summary <- NULL
for (each_tumor in tumor_type) {
  current <- match(each_tumor, tumor_type,nomatch = 0)
  print(paste("Processing the ", current, "tumors with total numbers of tumor", length(tumor_type),sep = " "))
  each_tumor_info_sum <- NULL
  each_tumor_info <- amp_delete[tumor_type == each_tumor, ]
  each_tumor_gene <- each_tumor_info$gene_mutation
  each_tumor_uniq_gene <- unique(each_tumor_gene)
  tumor_type_sum <- character(length(each_tumor_gene))
  gene_mut_sum <- character(length(each_tumor_gene))
  p_value_sum <- numeric(length(each_tumor_gene))
  numbersamples_withgene_sum <-numeric(length(each_tumor_gene))
  numbersamples_withoutgene_sum <- numeric(length(each_tumor_gene))
  total_others_sum <- numeric(length(each_tumor_gene))
  total_others_without_sum <- numeric(length(each_tumor_gene))
  each_tumor_total_sample <-   length(unique(each_tumor_info$sample_names))
  for (index in 1:length(each_tumor_uniq_gene)) {
    #current <- match(each_gene, each_tumor_gene,nomatch = 0)
    print(paste("Processing the ", index, "gene mutation, with total gene mutations", length(each_tumor_uniq_gene),sep = " "))
    each_gene <- each_tumor_uniq_gene[index]
    total_others <- 0
    total_others_without <- 0
    sample_loc <- which(each_gene == each_tumor_gene)
       #length(unique(each_tumor_info[gene_mutation == each_gene, ]$sample_names))
    numbersamples_withgene <- length(unique(each_tumor_info$sample_names[sample_loc]))
    #length(which(each_gene == each_tumor_gene))
      #length(unique(each_tumor_info[gene_mutation == each_gene, ]$sample_names))
    
    numbersamples_withoutgene <- each_tumor_total_sample- numbersamples_withgene
    
    outside_target_gene <- each_tumor_uniq_gene[!each_tumor_uniq_gene %in% each_gene]
    
    for (other_gene in outside_target_gene) {
 
      #nubmersamples_others_variants <- 
       other_sample_loc <- which(other_gene == each_tumor_gene)
       nubmersamples_others_variants <- length(unique(each_tumor_info$sample_names[other_sample_loc]))
      #ength(unique(each_tumor_info[gene_mutation == other_gene , .(sample_names)]$sample_names))
      
      #print(c(other_gene, nubmersamples_others_variants,old_nubmersamples_others_variants))
      numersamples_others_withoutvariants <- each_tumor_total_sample - nubmersamples_others_variants
      
      total_others <- total_others + nubmersamples_others_variants
      total_others_without <- total_others_without + numersamples_others_withoutvariants
    }
    testor <-rbind(c(numbersamples_withgene,numbersamples_withoutgene),
                   c(total_others, total_others_without))
    
    each_gene_p_value <- fisher.test(testor, alternative = "greater")$p.value
    
    # info <- data.table(tumor_type = each_tumor,
    #                    gene_mut = each_gene,
    #                    p_value = each_gene_p_value,
    #                    numbersamples_withgene = numbersamples_withgene,
    #                    numbersamples_withoutgene = numbersamples_withoutgene,
    #                    total_others = total_others,
    #                    total_others_without = total_others_without)
    # 
    # each_tumor_info_sum <- rbindlist(list(each_tumor_info_sum, info))
    tumor_type_sum[index] <- each_tumor
    gene_mut_sum[index] <- each_gene
    p_value_sum[index] <- each_gene_p_value
    numbersamples_withgene_sum[index] <- numbersamples_withgene
    numbersamples_withoutgene_sum[index] <- numbersamples_withoutgene
    total_others_without_sum[index] <- total_others_without
    total_others_sum[index] <- total_others
    
    }
  each_tumor_info_sum <- data.table(
    tumor_type = tumor_type_sum,
    gene_mut = gene_mut_sum,
    p_value = p_value_sum,
    numbersamples_withgene = numbersamples_withgene_sum,
    numbersamples_withoutgene = numbersamples_withoutgene_sum,
    total_others_withoutgene = total_others_without_sum,
    total_others = total_others_sum)
  
    #each_tumor_info_sum <- each_tumor_info_sum[order(p_value)]
    #each_tumor_info_sum$BH_pvalue = p.adjust(each_tumor_info_sum$p_value, method = "BH")
    
    total_gene_summary <- rbindlist(list(total_gene_summary, each_tumor_info_sum))
}
