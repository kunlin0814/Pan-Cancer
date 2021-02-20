library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")

base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_18.txt") 
exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

mutect_after_vaf <- fread(paste(base_dir,"NonSyn_Burair_filtering3_WithBreeds_Subtypes_QCpass_mutect_after_vaf_02_11.txt",
                                sep =seperator))
## indel
mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude,]
indel_file <- fread(paste(base_dir,"passQC_pan-tumor-total_indel_info_0214.txt",sep =seperator))
indel_file <- indel_file[gene_name!="-" & status=="frameshift" & ! sample_names %in% exclude,]
setcolorder(indel_file,c("sample_names","gene_name","emsembl_id","status"))
indel_file <- indel_file[,emsembl_id:=NULL]

Subtype <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table )
indel_file$Subtype <- Subtype
SNV <- unique(mutect_after_vaf[,c("sample_names","gene_name","status","Subtype"), with =F])
total_mut <- rbindlist(list(SNV,indel_file))
total_mut <- total_mut[!sample_names %in% exclude,]
## identify candidate genes (gene is 10 samples of all tumor and 5 samples within each tumor)

all_tumor_cut <- 10
signle_tumor_cut <- 5

all_gene <- unique(total_mut$gene_name)
all_tumor_type <- unique(total_mut$Subtype)
# might only have one sample has that
total_tumor_gene_sum <- NULL
for( index in 1:length(all_tumor_type)){
  print(paste("processing the",index,"tumor, with total tumors",length(all_tumor_type),sep = " " ))
  each_tumor <- all_tumor_type[index]
  each_tumor_gene_info <- unique(total_mut[Subtype == each_tumor, .(gene_name)][["gene_name"]])
  each_tumor_gene_candidate <- NULL
  target_gene_sample_number <- NULL
  target_gene_total_sample_number <- NULL
  for (gene_index in 1:length(each_tumor_gene_info)){
    each_gene <- each_tumor_gene_info[gene_index]
    #each_gene <- "PIK3CA"
    total_sample <- unique(total_mut[gene_name==each_gene, .(sample_names,Subtype)])
    total_sample_sum <- as.data.table(table(total_sample$Subtype))
    total_sample_number <- sum(total_sample_sum$N)
    each_tumor_sample_sum <- total_sample_sum[which(total_sample_sum$V1==each_tumor)]$N
    if (each_tumor_sample_sum > signle_tumor_cut && total_sample_number> all_tumor_cut){
      target_gene_sample_number <- c(target_gene_sample_number,each_tumor_sample_sum)
      target_gene_total_sample_number <- c(target_gene_total_sample_number,total_sample_number)
      each_tumor_gene_candidate <- c(each_tumor_gene_candidate,each_gene)
      }
    }
    each_tumor_sum <-  data.table(Subtype = each_tumor,
                                  gene_name = each_tumor_gene_candidate,
                                  target_gene_sample_number = target_gene_sample_number,
                                  target_gene_total_sample_number = target_gene_total_sample_number) 
    
    total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_tumor_sum),fill = T)
}
fwrite(total_tumor_gene_sum, file = "C:/Users/abc73/Desktop/candidate_gene_associated_TMB.txt",
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

## excldue NA and UCL
tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")


total_tumor_gene_sum <- na.omit(total_tumor_gene_sum[Subtype!="UCL"])
candidate_gene <- unique(total_tumor_gene_sum$gene_name)
## append TMB info
TMB_info <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
colnames(TMB_info)
total_mut$tmb <- match_vector_table(total_mut$sample_names, "combine_snv_indel_tmb", TMB_info, string_value = F)

# each tumor type seperate tmbl and tmbh
total_gene_sum <- NULL
for (each_gene in candidate_gene){
  tmb_l_group <- total_mut[Subtype %in%tmb_l, ]
  tmb_l_group_total_samples <- unique(tmb_l_group$sample_names)
  tmb_l_gene_mut_samples <- unique(tmb_l_group[gene_name==each_gene,.(sample_names)][["sample_names"]])
  tmb_l_gene_no_mut_samples <- setdiff(tmb_l_group_total_samples,tmb_l_gene_mut_samples)
  tmb_l_gene_mut_tmb <- unique(total_mut[sample_names %in% tmb_l_gene_mut_samples,.(sample_names,tmb)])[["tmb"]]
  tmb_l_gene_no_mut_tmb <- unique(total_mut[sample_names %in% tmb_l_gene_no_mut_samples,.(sample_names,tmb)])[["tmb"]]
  median_tmb_l_gene_mut_tmb <- median(tmb_l_gene_mut_tmb)
  median_tmb_l_gene_no_mut_tmb <- median(tmb_l_gene_no_mut_tmb)
  tmb_l_test <- wilcox.test(tmb_l_gene_mut_tmb,tmb_l_gene_no_mut_tmb)
  tmb_l_pvalue <- tmb_l_test$p.value
  
  tmb_h_group <- total_mut[Subtype %in%tmb_h, ]
  tmb_h_group_total_samples <- unique(tmb_h_group$sample_names)
  tmb_h_gene_mut_samples <- unique(tmb_h_group[gene_name==each_gene,.(sample_names)][["sample_names"]])
  tmb_h_gene_no_mut_samples <- setdiff(tmb_h_group_total_samples,tmb_h_gene_mut_samples)
  tmb_h_gene_mut_tmb <- unique(total_mut[sample_names %in% tmb_h_gene_mut_samples,.(sample_names,tmb)])[["tmb"]]
  tmb_h_gene_no_mut_tmb <- unique(total_mut[sample_names %in% tmb_h_gene_no_mut_samples,.(sample_names,tmb)])[["tmb"]]
  tmb_h_test <- wilcox.test(tmb_h_gene_mut_tmb,tmb_h_gene_no_mut_tmb)
  tmb_h_pvalue <- tmb_h_test$p.value
  median_tmb_h_gene_mut_tmb <- median(tmb_h_gene_mut_tmb)
  median_tmb_h_gene_no_mut_tmb <- median(tmb_h_gene_no_mut_tmb)
  
  each_gene_sum <- list(gene = each_gene,
                        median_tmb_l_gene_mut_tmb=median_tmb_l_gene_mut_tmb,
                        median_tmb_l_gene_no_mut_tmb=median_tmb_l_gene_no_mut_tmb,
                        median_tmb_h_gene_mut_tmb= median_tmb_h_gene_mut_tmb,
                        median_tmb_h_gene_no_mut_tmb=median_tmb_h_gene_no_mut_tmb,
                        tmb_l_pvalue=tmb_l_pvalue,
                        tmb_h_pvalue=tmb_h_pvalue)
  
  
  total_gene_sum <- rbindlist(list(total_gene_sum,each_gene_sum))
}
fwrite(total_tumor_gene_sum, file = "C:/Users/abc73/Desktop/p_value_candidate_gene_associated_TMB.txt",
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")
# 
# 
# tmb_h_group <- total_mut[Subtype %in%tmn_h, ]
# tmb_l_gene_mut_group <- 
# tmb_l_gene_no_mut_group <- 
# 
# 
