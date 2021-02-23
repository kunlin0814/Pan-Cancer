library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")

base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

output_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/Gene_association_TMB/"

pathway <- fread(paste(base_dir,"all_pathway.txt",sep = seperator), na.strings = "")
target_pathway_gene <- fread(paste(base_dir,"target_pathway_total_genes.txt",sep = seperator), na.strings = "")

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 
exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

mutect_after_vaf <- fread(paste(base_dir,"NonSyn_Burair_filtering3_WithBreeds_Subtypes_QCpass_mutect_after_vaf_02_11.txt",
                                sep =seperator))

s1_data <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
exclude <- c(exclude, s1_high_sample)

## indel
mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude,]
indel_file <- fread(paste(base_dir,"passQC_pan-tumor-total_indel_info_0214.txt",sep =seperator))
indel_file <- indel_file[gene_name!="-" & status=="frameshift" & ! sample_names %in% exclude,]
setcolorder(indel_file,c("sample_names","gene_name","emsembl_id","status"))
indel_file <- indel_file[,emsembl_id:=NULL]

Subtype <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table )
indel_file$Subtype <- Subtype
SNV <- unique(mutect_after_vaf[,c("sample_names","gene_name","status","Subtype"), with =F])


amp_delete <- fread(paste(base_dir,"CNV_Taifang_total_amp_delete_no_pseudo_subtype_02_18.txt",sep = seperator),
                    header = T,na.strings = "")

amp_delete <- amp_delete[!sample_names %in% exclude]
amp_delete <- amp_delete[!grepl("ENSCAFG",amp_delete[,.(gene_name)]$gene_name,ignore.case = T)]

amp_delete_data <- amp_delete[gene_name %in% target_pathway_gene$target_pathway_gene,.(sample_names,gene_name,mut_type,subtype)]
colnames(amp_delete_data) <- c("sample_names","gene_name","status","Subtype")
total_mut <- rbindlist(list(SNV,indel_file,amp_delete_data))
total_mut <- total_mut[!sample_names %in% exclude & Subtype !="UCL",]

## identify candidate genes (gene is 10 samples of all tumor and 5 samples within each tumor)
## need to normalize with median for each tumor
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
    # each_gene <- "PIK3CA"
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
fwrite(total_tumor_gene_sum, file = paste(output_dir,"02_19","include_amp_candidate_gene_associated_TMB_02_22.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

## Use candidate gene to compare tmb
base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_18.txt") 
exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

s1_data <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
exclude <- c(exclude, s1_high_sample)

## excldue NA and UCL
tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")

total_tumor_gene_sum <- fread(paste(output_dir,"02_19","include_amp_candidate_gene_associated_TMB_02_22.txt",sep = seperator))

total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)
candidate_gene <- unique(total_tumor_gene_sum$gene_name)
## append TMB info
TMB_info <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
colnames(TMB_info)
total_mut$tmb <- match_vector_table(total_mut$sample_names, "combine_snv_indel_tmb", TMB_info, string_value = F)

## Normalize TMB with regards to each tumor median (x-u/std)
total_tumor_type <- unique(total_mut$Subtype)
total_tumor_normalize <- NULL

for (each_tumor in total_tumor_type){
  each_tumor_info <- total_mut[Subtype==each_tumor,]
  each_median <- median(total_mut[Subtype==each_tumor, .(tmb)][['tmb']])
  # each_sd <- sd(total_mut[Subtype==each_tumor, .(tmb)][['tmb']]) 
  each_tumor_info <- each_tumor_info[, normalizetmb:= log10((tmb/each_median)+0.1)]
  total_tumor_normalize <- rbindlist(list(total_tumor_normalize,each_tumor_info))
}

total_mut <- total_tumor_normalize


# each_median <- median(total_mut[Subtype=="MT", .(tmb)][['tmb']])
# each_tumor_info <- total_mut[Subtype=="MT", .(tmb)]
# # each_sd <- sd(total_mut[Subtype==each_tumor, .(tmb)][['tmb']]) 
# each_tumor_info[, normalizetmb:= log10((tmb/each_median)+0.1)]
# 
# log10(each_tumor_info$tmb/each_median+0.1)
## Normalize end 


# each tumor type seperate tmbl and tmbh
total_gene_sum <- NULL
for (each_gene in candidate_gene){
  #print(each_gene)
  #each_gene <- "TRAV19"
  tmb_l_group <- total_mut[Subtype %in%tmb_l, ]
  tmb_l_group_total_samples <- unique(tmb_l_group$sample_names)
  tmb_l_gene_mut_samples <- unique(tmb_l_group[gene_name==each_gene,.(sample_names)][["sample_names"]])
  tmb_l_gene_no_mut_samples <- setdiff(tmb_l_group_total_samples,tmb_l_gene_mut_samples)
  tmb_l_gene_mut_tmb <- unique(total_mut[sample_names %in% tmb_l_gene_mut_samples,.(sample_names,normalizetmb)])[["normalizetmb"]]
  tmb_l_gene_no_mut_tmb <- unique(total_mut[sample_names %in% tmb_l_gene_no_mut_samples,.(sample_names,normalizetmb)])[["normalizetmb"]]
  
  
  if (length(tmb_l_gene_mut_samples)!=0 && length(tmb_l_gene_no_mut_samples)!=0){
    tmb_l_test <- wilcox.test(tmb_l_gene_mut_tmb,tmb_l_gene_no_mut_tmb)
    tmb_l_pvalue <- tmb_l_test$p.value
    median_tmb_l_gene_mut_tmb <- median(tmb_l_gene_mut_tmb)
    median_tmb_l_gene_no_mut_tmb <- median(tmb_l_gene_no_mut_tmb)
  }
  else{
    tmb_l_pvalue <- "No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
    median_tmb_l_gene_mut_tmb <- "No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
    median_tmb_l_gene_no_mut_tmb <-"No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
  }
  
  tmb_h_group <- total_mut[Subtype %in%tmb_h, ]
  tmb_h_group_total_samples <- unique(tmb_h_group$sample_names)
  tmb_h_gene_mut_samples <- unique(tmb_h_group[gene_name==each_gene,.(sample_names)][["sample_names"]])
  tmb_h_gene_no_mut_samples <- setdiff(tmb_h_group_total_samples,tmb_h_gene_mut_samples)
  tmb_h_gene_mut_tmb <- unique(total_mut[sample_names %in% tmb_h_gene_mut_samples,.(sample_names,normalizetmb)])[["normalizetmb"]]
  tmb_h_gene_no_mut_tmb <- unique(total_mut[sample_names %in% tmb_h_gene_no_mut_samples,.(sample_names,normalizetmb)])[["normalizetmb"]]
  
  if (length(tmb_h_gene_mut_samples)!=0 && length(tmb_h_gene_no_mut_samples)!=0){
    
    tmb_h_test <- wilcox.test(tmb_h_gene_mut_tmb,tmb_h_gene_no_mut_tmb)
    tmb_h_pvalue <- tmb_h_test$p.value
    median_tmb_h_gene_mut_tmb <- median(tmb_h_gene_mut_tmb)
    median_tmb_h_gene_no_mut_tmb <- median(tmb_h_gene_no_mut_tmb)
  }
  
  else{
    tmb_h_pvalue <-"No tmb_h_gene_mut_samples or tmb_h_gene_no_mut_samples"
    median_tmb_h_gene_mut_tmb <- "No tmb_h_gene_mut_samples or tmb_h_gene_no_mut_samples"
    median_tmb_h_gene_no_mut_tmb <-"No tmb_h_gene_mut_samples or tmb_h_gene_no_mut_samples"
  }
  
  # if (length(tmb_h_gene_mut_samples)!=0 && length(tmb_h_gene_no_mut_samples)!=0 
  #     && length(tmb_l_gene_mut_samples)!=0 && length(tmb_l_gene_no_mut_samples)!=0){
  each_gene_sum <- list(gene = each_gene,
                        median_tmb_l_gene_mut_tmb=median_tmb_l_gene_mut_tmb,
                        median_tmb_l_gene_no_mut_tmb=median_tmb_l_gene_no_mut_tmb,
                        median_tmb_h_gene_mut_tmb= median_tmb_h_gene_mut_tmb,
                        median_tmb_h_gene_no_mut_tmb=median_tmb_h_gene_no_mut_tmb,
                        tmb_l_pvalue=tmb_l_pvalue,
                        tmb_h_pvalue=tmb_h_pvalue)
  
  
  total_gene_sum <- rbindlist(list(total_gene_sum,each_gene_sum))
}

fwrite(total_gene_sum, file = paste(output_dir,"02_19","include_amp_normalize_p_value_candidate_gene_associated_TMB_02_22.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")
# 
# 
# tmb_h_group <- total_mut[Subtype %in%tmn_h, ]
# tmb_l_gene_mut_group <- 
# tmb_l_gene_no_mut_group <- 
# 
# 
