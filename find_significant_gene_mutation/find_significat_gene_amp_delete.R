library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/Burair_filtering3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
seperator <- "/"

#retro_gene_list <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Retro_gene_finding/RetroGeneList/new_retro_gene_list_CanFam3.1.99gtf.txt",
#                         header = F)
whole_wes_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt") 
#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt")       

exclude <- unique(unlist(whole_wes_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

amp_delete <- fread(paste(base_dir,"total_samples_amp_delete.txt",sep = seperator),header = F)
colnames(amp_delete) <- c("sample_names","gene_name","mutation_type","CNA","tumor_type")
amp_delete <- amp_delete[,gene_mutation:=paste(gene_name,mutation_type,sep = "_")]


tumor_type <- sort(unique(amp_delete$tumor_type))
tumor_type <- "GLM"

total_gene_summary <- NULL
for (each_tumor in tumor_type) {
  current <- match(each_tumor, tumor_type,nomatch = 0)
  print(paste("Processing the ", current, "tumors with total numbers of tumor", length(tumor_type),sep = " "))
  each_tumor_info_sum <- NULL
  each_tumor_info <- amp_delete[tumor_type == each_tumor, ]
  each_tumor_gene <- unique(each_tumor_info$gene_mutation)
  
  for (each_gene in each_tumor_gene) {
    current <- match(each_gene, each_tumor_gene,nomatch = 0)
    print(paste("Processing the ", current, "gene mutation, with total gene mutations", length(each_tumor_gene),sep = " "))
    
    total_others <- 0
    total_others_without <- 0

    numbersamples_withgene <- length(unique(each_tumor_info[gene_mutation == each_gene, ][["sample_names"]]))
    
    numbersamples_withoutgene <- length(unique(each_tumor_info$sample_names)) - numbersamples_withgene
    
    outside_target_gene <-each_tumor_gene[!each_tumor_gene %in% each_gene]
    
    for (other_gene in outside_target_gene) {
      nubmersamples_others_variants <-
        length(unique(each_tumor_info[gene_mutation == other_gene , .(sample_names)][['sample_names']]))
      
      numersamples_others_withoutvariants <-
        length(unique(each_tumor_info$sample_names)) - nubmersamples_others_variants
      
      total_others <- total_others + nubmersamples_others_variants
      total_others_without <- total_others_without + numersamples_others_withoutvariants
    }
    testor <-rbind(c(numbersamples_withgene,numersamples_withoutgene),
                   c(total_others, total_others_without))
    
    each_gene_p_value <- fisher.test(testor, alternative = "greater")$p.value
    
    info <- data.table(tumor_type = each_tumor,
                       gene_mut = each_gene,
                       p_value = each_gene_p_value,
                       numbersamples_withgene = numbersamples_withgene,
                       numbersamples_withoutgene = numbersamples_withoutgene,
                       total_others = total_others,
                       total_others_without = total_others_without)
    
    each_tumor_info_sum <- rbindlist(list(each_tumor_info_sum, info))
    }
  
    each_tumor_info_sum <- each_tumor_info_sum[order(p_value)]
    each_tumor_info_sum$BH_pvalue = p.adjust(each_tumor_info_sum$p_value, method = "BH")
    
    total_gene_summary <- rbindlist(list(total_gene_summary, each_tumor_info_sum))
}
