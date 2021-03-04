library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")


## Gene assocaited TMB only exam somatic mutation
base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

output_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/Gene_association_TMB"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/Gene_association_TMB/"

whole_wes_clean_breed_table <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt")
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 
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


## exclude s1 high and UCL samples
s1_data <- fread("/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
exclude <- c(exclude, s1_high_sample)

total_mut <- total_mut[!sample_names %in% exclude & Subtype !="UCL",]
## append TMB 
TMB_info <- fread("/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
colnames(TMB_info)
total_mut$tmb <- match_vector_table(total_mut$sample_names, "combine_snv_indel_tmb", TMB_info, string_value = F)

## Normalize TMB with regards to each tumor median (in pan-tumor analysis)
## now decide not normalize 2/25
total_tumor_type <- unique(total_mut$Subtype)
total_tumor_normalize <- NULL

for (each_tumor in total_tumor_type){
  each_tumor_info <- total_mut[Subtype==each_tumor,]
  each_median <- median(total_mut[Subtype==each_tumor, .(tmb)][['tmb']])
  # each_sd <- sd(total_mut[Subtype==each_tumor, .(tmb)][['tmb']])
  each_tumor_info <- each_tumor_info[, normalizetmb:= (tmb/each_median)+0.1]
  #each_tumor_info <- each_tumor_info[, normalizetmb:= tmb]
  total_tumor_normalize <- rbindlist(list(total_tumor_normalize,each_tumor_info))
}

total_mut <- total_tumor_normalize

## Normalize end 

## identify candidate genes (5 samples within each tumor)
## need to normalize with median for each tumor in pan tumor analysis

## Pan_tumor cross tumor types
compare_gene <- "TP53"
signle_tumor_cut <- 5
total_tumor_gene_sum <- NULL
pan_tumor_uniq_gene <- unique(total_mut$gene_name)
total_sample_number <- length(unique(total_mut$sample_names))
for (gene_index in 1:length(pan_tumor_uniq_gene)){
  each_gene_summary <- list()
  each_gene <- pan_tumor_uniq_gene[gene_index]
  #each_gene <- "PIK3CA"
  each_gene_total_sample <- unique(total_mut[gene_name==each_gene, .(sample_names)])
  each_gene_total_sample_number <- nrow(each_gene_total_sample)
  if (each_gene_total_sample_number >= signle_tumor_cut ){
    candidate_gene <- each_gene
    gene_mut_tmb <- unique(total_mut[gene_name==candidate_gene, .(sample_names,normalizetmb)])[['normalizetmb']]
    gene_mut_sample <- unique(total_mut[gene_name==candidate_gene, .(sample_names)])[['sample_names']]
    
    gene_no_mut_tmb <-  unique(total_mut[!sample_names %in% gene_mut_sample, .(sample_names,normalizetmb)])[['normalizetmb']]
    tmb_test <- wilcox.test(gene_mut_tmb,gene_no_mut_tmb)
    p_value <- tmb_test$p.value
    fold_change <- median(gene_mut_tmb)/median(gene_no_mut_tmb)
    
    compare_gene_sample <- unique(total_mut[gene_name==compare_gene,][['sample_names']])
    candidate_gene_sample<- unique(total_mut[gene_name==candidate_gene][['sample_names']])
    alt_alt <- length(intersect(compare_gene_sample,candidate_gene_sample)) #TP53 and another gene sample
    alt_no_alt <- length(setdiff(compare_gene_sample,candidate_gene_sample)) # TP53 but not other gene mut sample
    no_alt_alt <- length(setdiff(candidate_gene_sample,compare_gene_sample))# other gene mut sample but not tp53 mut sample
    no_alt_no_alt <- total_sample_number-alt_alt-alt_no_alt-no_alt_alt
    tp53_fisher_test <- fisher.test(rbind(c(alt_alt, alt_no_alt), c(no_alt_alt, no_alt_no_alt))); 
    tp53_fisher_p <- tp53_fisher_test[["p.value"]];
    tp53_relation_type <- ifelse(tp53_fisher_test[["estimate"]] > 1, "Inclusive", "Exclusive");
    tp53_relation_sign <- ifelse(tp53_fisher_p <= 0.05, "Yes", "No");
    each_gene_summary <- list(Gene = candidate_gene,
                              Mutated_samples = each_gene_total_sample_number,
                              P_value=p_value,
                              Fold_change= fold_change,
                              TP53_mutual_P_value=tp53_fisher_p,
                              TP53_Incl_Excl=tp53_relation_type,
                              TP53_mutual_significant = tp53_relation_sign)
    
    
  }
  total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_gene_summary))
}
total_tumor_gene_sum <- setDT(total_tumor_gene_sum)
total_tumor_gene_sum <- total_tumor_gene_sum[order(P_value)]
total_tumor_gene_sum$BH_P_value <-  p.adjust(total_tumor_gene_sum$P_value, method = "BH")

fwrite(total_tumor_gene_sum, file = paste(output_dir,"03_02","pan_tumor_gene_assication_TMB.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")
## cross tumor end

### examine each tumor type
total_tumor_type <- unique(total_mut$Subtype)
compare_gene <- "TP53"
signle_tumor_cut <- 5
total_tumor_gene_sum <- NULL
total_sample_number <- length(unique(total_mut$sample_names))
for (tumor_index in 1:length(total_tumor_type)){
  each_tumor <- total_tumor_type[tumor_index]
  each_tumor_info <- total_mut[Subtype==each_tumor,]
  each_tumor_uniq_gene <- unique(each_tumor_info$gene_name)
  each_tumor_total_gene_summary <- NULL
  for (gene_index in 1:length(each_tumor_uniq_gene)){
  each_tumor_each_gene_summary <- NULL
  each_gene <- each_tumor_uniq_gene[gene_index]
  #print(each_gene)
  #each_gene <- "PIK3CA"
  each_tumor_gene_total_sample <- unique(each_tumor_info[gene_name==each_gene, .(sample_names)])
  each_tumor_gene_total_sample_number <- nrow(each_tumor_gene_total_sample)
  if (each_tumor_gene_total_sample_number >= signle_tumor_cut ){
    candidate_gene <- each_gene
    gene_mut_tmb <- unique(each_tumor_info[gene_name==candidate_gene, .(sample_names,tmb)])[['tmb']]
    gene_mut_sample <- unique(each_tumor_info[gene_name==candidate_gene, .(sample_names)])[['sample_names']]
    
    gene_no_mut_tmb <-  unique(each_tumor_info[!sample_names %in% gene_mut_sample, .(sample_names,tmb)])[['tmb']]
    tmb_test <- wilcox.test(gene_mut_tmb,gene_no_mut_tmb)
    p_value <- tmb_test$p.value
    fold_change <- median(gene_mut_tmb)/median(gene_no_mut_tmb)
    
    compare_gene_sample <- unique(each_tumor_info[gene_name==compare_gene,][['sample_names']])
    candidate_gene_sample<- unique(each_tumor_info[gene_name==candidate_gene][['sample_names']])
    alt_alt <- length(intersect(compare_gene_sample,candidate_gene_sample)) #TP53 and another gene sample
    alt_no_alt <- length(setdiff(compare_gene_sample,candidate_gene_sample)) # TP53 but not other gene mut sample
    no_alt_alt <- length(setdiff(candidate_gene_sample,compare_gene_sample))# other gene mut sample but not tp53 mut sample
    no_alt_no_alt <- total_sample_number-alt_alt-alt_no_alt-no_alt_alt
    tp53_fisher_test <- fisher.test(rbind(c(alt_alt, alt_no_alt), c(no_alt_alt, no_alt_no_alt))); 
    tp53_fisher_p <- tp53_fisher_test[["p.value"]];
    tp53_relation_type <- ifelse(tp53_fisher_test[["estimate"]] > 1, "Inclusive", "Exclusive");
    tp53_relation_sign <- ifelse(tp53_fisher_p <= 0.05, "Yes", "No");
    each_tumor_each_gene_summary <- list(Tumor_type= each_tumor,
                              Gene = candidate_gene,
                              Mutated_samples = each_tumor_gene_total_sample_number,
                              P_value=p_value,
                              Fold_change= fold_change,
                              TP53_mutual_P_value=tp53_fisher_p,
                              TP53_Incl_Excl=tp53_relation_type,
                              TP53_mutual_significant = tp53_relation_sign)
    
    
  }
  each_tumor_total_gene_summary <- rbindlist(list(each_tumor_total_gene_summary,each_tumor_each_gene_summary))  
  
  }
  each_tumor_total_gene_summary <- setDT(each_tumor_total_gene_summary)
  each_tumor_total_gene_summary <- each_tumor_total_gene_summary[order(P_value)]
  each_tumor_total_gene_summary$BH_P_value <-  p.adjust(each_tumor_total_gene_summary$P_value, method = "BH")
  total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_tumor_total_gene_summary))
  }
 
fwrite(total_tumor_gene_sum, file = paste(output_dir,"03_02","all_tumor_type_gene_assication_TMB.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

## TMB-l and TMB-h gene associated tmb

## Gene assocaited TMB only exam somatic mutation
base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

output_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/Gene_association_TMB"
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/Gene_association_TMB/"

whole_wes_clean_breed_table <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt")
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 
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

## exclude s1 high and UCL samples
s1_data <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
exclude <- c(exclude, s1_high_sample)

total_mut <- total_mut[!sample_names %in% exclude & Subtype !="UCL",]
## append TMB 
TMB_info <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
colnames(TMB_info)
total_mut$tmb <- match_vector_table(total_mut$sample_names, "combine_snv_indel_tmb", TMB_info, string_value = F)

all_tumor_cut <- 0
signle_tumor_cut <- 5
tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")

### TMB-l candidate genes
# might only have one sample has that
total_tumor_gene_sum <- NULL
for( index in 1:length(tmb_l)){
  print(paste("processing the",index,"tumor, with total tumors",length(tmb_l),sep = " " ))
  each_tumor <- tmb_l[index]
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
    if (each_tumor_sample_sum >= signle_tumor_cut && total_sample_number>= all_tumor_cut){
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
total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)
fwrite(total_tumor_gene_sum, file = paste(output_dir,"03_02","tmb_l_Not_include_amp_candidate_gene_associated_TMB_03_02.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

### TMB-h candidate genes

# might only have one sample has that
total_tumor_gene_sum <- NULL
for( index in 1:length(tmb_h)){
  print(paste("processing the",index,"tumor, with total tumors",length(tmb_h),sep = " " ))
  each_tumor <- tmb_h[index]
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
    if (each_tumor_sample_sum >= signle_tumor_cut && total_sample_number>= all_tumor_cut){
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
total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)
fwrite(total_tumor_gene_sum, file = paste(output_dir,"03_02","tmb_h_Not_include_amp_candidate_gene_associated_TMB_03_02.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

## Use candidate gene to compare tmb
tmb_l_total_tumor_gene_sum <- fread(paste(output_dir,"03_02","tmb_l_Not_include_amp_candidate_gene_associated_TMB_03_02.txt",sep = seperator))
tmb_l_total_tumor_gene_sum <- na.omit(tmb_l_total_tumor_gene_sum)
## excldue NA and UCL
tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")
compare_gene <- "TP53"
total_mut <- total_mut[!sample_names %in% exclude & Subtype !="UCL",]
## append TMB 
TMB_info <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")

tmb_l_candidate_genes <- unique(tmb_l_total_tumor_gene_sum$gene_name)
## append TMB info
colnames(TMB_info)
total_mut$tmb <- match_vector_table(total_mut$sample_names, "combine_snv_indel_tmb", TMB_info, string_value = F)

## Normalize TMB with regards to each tumor median
## now decide not normalize 2/25
# each tumor type seperate tmbl and tmbh
# tmb-l
tmb_l_total_gene_sum <- NULL
for (each_gene in tmb_l_candidate_genes){
  print(each_gene)
  tmb_l_group <- total_mut[Subtype %in%tmb_l, ]
  tmb_l_group_total_samples <- unique(tmb_l_group$sample_names)
  tmb_l_gene_mut_samples <- unique(tmb_l_group[gene_name==each_gene,.(sample_names)][["sample_names"]])
  tmb_l_gene_no_mut_samples <- unique(setdiff(tmb_l_group_total_samples,tmb_l_gene_mut_samples))
  tmb_l_gene_mut_tmb <- unique(total_mut[sample_names %in% tmb_l_gene_mut_samples,.(sample_names,tmb)])[["tmb"]]
  tmb_l_gene_no_mut_tmb <- unique(total_mut[sample_names %in% tmb_l_gene_no_mut_samples,.(sample_names,tmb)])[["tmb"]]
  each_gene_summary <- list()
  if (length(tmb_l_gene_mut_samples)!=0 || length(tmb_l_gene_no_mut_samples)!=0){
    candidate_gene <- each_gene
    tmb_l_group_total_samples_number <- length(tmb_l_group_total_samples)
    median_tmb_l_gene_mut_tmb <- median(tmb_l_gene_mut_tmb)
    median_tmb_l_gene_no_mut_tmb <- median(tmb_l_gene_no_mut_tmb)
    tmb_l_fold_change <- median(median_tmb_l_gene_mut_tmb)/median(median_tmb_l_gene_no_mut_tmb)
    tmb_l_test <- wilcox.test(tmb_l_gene_mut_tmb,tmb_l_gene_no_mut_tmb)
    tmb_l_pvalue <- tmb_l_test$p.value
    
    tmb_l_compare_gene_sample <- unique(tmb_l_group[gene_name==compare_gene,][['sample_names']])
    tmb_l_candidate_gene_sample<- unique(tmb_l_group[gene_name==candidate_gene][['sample_names']])
    tmb_l_alt_alt <- length(intersect(tmb_l_compare_gene_sample,tmb_l_candidate_gene_sample)) #TP53 and another gene sample
    tmb_l_alt_no_alt <- length(setdiff(tmb_l_compare_gene_sample,tmb_l_candidate_gene_sample)) # TP53 but not other gene mut sample
    tmb_l_no_alt_alt <- length(setdiff(tmb_l_candidate_gene_sample,tmb_l_compare_gene_sample))# other gene mut sample but not tp53 mut sample
    tmb_l_no_alt_no_alt <- tmb_l_group_total_samples_number-tmb_l_alt_alt-tmb_l_alt_no_alt-tmb_l_no_alt_alt
    tmb_l_tp53_fisher_test <- fisher.test(rbind(c(tmb_l_alt_alt, tmb_l_alt_no_alt), c(tmb_l_no_alt_alt, tmb_l_no_alt_no_alt))); 
    tmb_l_tp53_fisher_p <- tmb_l_tp53_fisher_test[["p.value"]];
    tmb_l_tp53_relation_type <- ifelse(tmb_l_tp53_fisher_test[["estimate"]] > 1, "Inclusive", "Exclusive");
    tmb_l_tp53_relation_sign <- ifelse(tmb_l_tp53_fisher_p <= 0.05, "Yes", "No");
    
    each_gene_summary <- list(
      Gene = candidate_gene,
      Mutated_samples = length(tmb_l_gene_mut_samples),
      P_value=tmb_l_pvalue,
      Fold_change= tmb_l_fold_change,
      TP53_mutual_P_value=tmb_l_tp53_fisher_p,
      TP53_Incl_Excl=tmb_l_tp53_relation_type,
      TP53_mutual_significant = tmb_l_tp53_relation_sign)
    
  }
  tmb_l_total_gene_sum <- rbindlist(list(tmb_l_total_gene_sum,each_gene_summary))
}

tmb_l_total_gene_sum <- setDT(tmb_l_total_gene_sum)
tmb_l_total_gene_sum <- tmb_l_total_gene_sum[order(P_value)]
tmb_l_total_gene_sum$tmb_l_BH_pvalue <-  p.adjust(tmb_l_total_gene_sum$P_value, method = "BH")
fwrite(tmb_l_total_gene_sum, file = paste(output_dir,"03_02","tmb_l_Not_include_amp_candidate_gene_associated_TMB_03_02_data_summary.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

## tmb-l end 


# tmb-h  
total_gene_sum <- NULL
for (each_gene in tmb_l_candidate_genes){
  print(each_gene)
  tmb_l_group <- total_mut[Subtype %in%tmb_l, ]
  tmb_l_group_total_samples <- unique(tmb_l_group$sample_names)
  tmb_l_gene_mut_samples <- unique(tmb_l_group[gene_name==each_gene,.(sample_names)][["sample_names"]])
  tmb_l_gene_no_mut_samples <- unique(setdiff(tmb_l_group_total_samples,tmb_l_gene_mut_samples))
  tmb_l_gene_mut_tmb <- unique(total_mut[sample_names %in% tmb_l_gene_mut_samples,.(sample_names,tmb)])[["tmb"]]
  tmb_l_gene_no_mut_tmb <- unique(total_mut[sample_names %in% tmb_l_gene_no_mut_samples,.(sample_names,tmb)])[["tmb"]]
  
  if (length(tmb_l_gene_mut_samples)!=0 || length(tmb_l_gene_no_mut_samples)!=0){
    each_gene <- candidate_gene
    tmb_l_group_total_samples_number <- length(tmb_l_group_total_samples)
    median_tmb_l_gene_mut_tmb <- median(tmb_l_gene_mut_tmb)
    median_tmb_l_gene_no_mut_tmb <- median(tmb_l_gene_no_mut_tmb)
    tmb_l_fold_change <- median(median_tmb_l_gene_mut_tmb)/median(median_tmb_l_gene_no_mut_tmb)
    tmb_l_test <- wilcox.test(tmb_l_gene_mut_tmb,tmb_l_gene_no_mut_tmb)
    tmb_l_pvalue <- tmb_l_test$p.value
    
    tmb_l_compare_gene_sample <- unique(tmb_l_group[gene_name==compare_gene,][['sample_names']])
    tmb_l_candidate_gene_sample<- unique(tmb_l_group[gene_name==candidate_gene][['sample_names']])
    tmb_l_alt_alt <- length(intersect(tmb_l_compare_gene_sample,tmb_l_candidate_gene_sample)) #TP53 and another gene sample
    tmb_l_alt_no_alt <- length(setdiff(tmb_l_compare_gene_sample,tmb_l_candidate_gene_sample)) # TP53 but not other gene mut sample
    tmb_l_no_alt_alt <- length(setdiff(tmb_l_candidate_gene_sample,tmb_l_compare_gene_sample))# other gene mut sample but not tp53 mut sample
    tmb_l_no_alt_no_alt <- tmb_l_group_total_samples_number-tmb_l_alt_alt-tmb_l_alt_no_alt-tmb_l_no_alt_alt
    tmb_l_tp53_fisher_test <- fisher.test(rbind(c(tmb_l_alt_alt, tmb_l_alt_no_alt), c(tmb_l_no_alt_alt, tmb_l_no_alt_no_alt))); 
    tmb_l_tp53_fisher_p <- tmb_l_tp53_fisher_test[["p.value"]];
    tmb_l_tp53_relation_type <- ifelse(tmb_l_tp53_fisher_test[["estimate"]] > 1, "Inclusive", "Exclusive");
    tmb_l_tp53_relation_sign <- ifelse(tmb_l_tp53_fisher_p <= 0.05, "Yes", "No");
    
  }
  # else{
  #   tmb_l_pvalue <- "No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
  #   median_tmb_l_gene_mut_tmb <- "No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
  #   median_tmb_l_gene_no_mut_tmb <-"No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
  # }
}  
  
  each_gene_sum <- list(Gene = candidate_gene,
                    
                        tmb_l_Mutated_samples=length(tmb_l_gene_mut_samples),
                        
                        tmb_h_gene_number_samples = length(tmb_h_gene_mut_samples),
                        tmb_h_gene_no_mut_samples = length(tmb_h_gene_no_mut_samples),
                        tmb_l_fold_change = median_tmb_l_gene_mut_tmb/median_tmb_l_gene_no_mut_tmb,
                        tmb_h_fold_change = median_tmb_h_gene_mut_tmb/median_tmb_h_gene_no_mut_tmb,
                        tmb_l_pvalue = tmb_l_pvalue,
                        tmb_h_pvalue = tmb_h_pvalue)
  
  each_tumor_each_gene_summary <- list(Tumor_type= each_tumor,
                                       Gene = candidate_gene,
                                       Mutated_samples = each_tumor_gene_total_sample_number,
                                       P_value=p_value,
                                       Fold_change= fold_change,
                                       TP53_mutual_P_value=tp53_fisher_p,
                                       TP53_Incl_Excl=tp53_relation_type,
                                       TP53_mutual_significant = tp53_relation_sign)
  
  total_gene_sum <- rbindlist(list(total_gene_sum,each_gene_sum))
  }
}
total_gene_sum <- setDT(total_gene_sum)
total_gene_sum <- total_gene_sum[order(tmb_l_pvalue)]
total_gene_sum$tmb_l_BH_pvalue <-  p.adjust(total_gene_sum$tmb_l_pvalue, method = "BH")
total_gene_sum <- total_gene_sum[order(tmb_h_pvalue)]
total_gene_sum$tmb_h_BH_pvalue <- p.adjust(total_gene_sum$tmb_h_pvalue,method = "BH")


fwrite(total_gene_sum, file = paste(output_dir,"03_02","Not_include_amp_not_normalize_p_value_candidate_gene_associated_TMB_03_02.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

### now seperate snv and amp
base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 
exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

amp_delete <- fread(paste(base_dir,"CNV_Taifang_total_amp_delete_no_pseudo_subtype_02_18.txt",sep = seperator),
                    header = T,na.strings = "")

amp_delete <- amp_delete[!sample_names %in% exclude]
amp_delete <- amp_delete[!grepl("ENSCAFG",amp_delete[,.(gene_name)]$gene_name,ignore.case = T)]
CNV <- unique(amp_delete[,c("sample_names","gene_name","mut_type","subtype"),with =F])
colnames(CNV)<- c("sample_names","gene_name","status","Subtype")
CNV <- CNV[!sample_names %in% exclude & Subtype!="UCL",]
#amp_delete_data <- amp_delete[gene_name %in% target_pathway_gene$target_pathway_gene,.(sample_names,gene_name,mut_type,subtype)]
total_cnv <- rbindlist(list(CNV))

all_tumor_cut <- 10
signle_tumor_cut <- 5

all_tumor_type <- unique(total_cnv$Subtype)
# might only have one sample has that
total_tumor_gene_sum <- NULL
for( index in 1:length(all_tumor_type)){
  print(paste("processing the",index,"tumor, with total tumors",length(all_tumor_type),sep = " " ))
  each_tumor <- all_tumor_type[index]
  each_tumor_gene_info <- unique(total_cnv[Subtype == each_tumor, .(gene_name)][["gene_name"]])
  each_tumor_gene_candidate <- NULL
  target_gene_sample_number <- NULL
  target_gene_total_sample_number <- NULL
  for (gene_index in 1:length(each_tumor_gene_info)){
    each_gene <- each_tumor_gene_info[gene_index]
    # each_gene <- "PIK3CA"
    total_sample <- unique(total_cnv[gene_name==each_gene, .(sample_names,Subtype)])
    total_sample_sum <- as.data.table(table(total_sample$Subtype))
    total_sample_number <- sum(total_sample_sum$N)
    each_tumor_sample_sum <- total_sample_sum[which(total_sample_sum$V1==each_tumor)]$N
    if (each_tumor_sample_sum >= signle_tumor_cut && total_sample_number>= all_tumor_cut){
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
total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)
fwrite(total_tumor_gene_sum, file = paste(output_dir,"03_02","amp_candidate_gene_associated_TMB_03_02.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

## Use candidate gene to compare tmb

## exclude s1 high
# s1_data <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
# s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
# exclude <- c(exclude, s1_high_sample)

## excldue NA and UCL
tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")

# total_tumor_gene_sum <- fread(paste(output_dir,"03_02","Not_include_amp_candidate_gene_associated_TMB_03_02.txt",sep = seperator))
# 
total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)

## append TMB info
TMB_info <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
colnames(TMB_info)
total_cnv$tmb <- match_vector_table(total_cnv$sample_names, "combine_snv_indel_tmb", TMB_info, string_value = F)
total_tumor_type <- unique(total_cnv$Subtype)

## Normalize TMB with regards to each tumor median
## now decide not normalize 2/25
# total_tumor_normalize <- NULL
# 
# for (each_tumor in total_tumor_type){
#   each_tumor_info <- total_cnv[Subtype==each_tumor,]
#   each_median <- median(total_cnv[Subtype==each_tumor, .(tmb)][['tmb']])
#   # each_sd <- sd(total_cnv[Subtype==each_tumor, .(tmb)][['tmb']]) 
#   each_tumor_info <- each_tumor_info[, normalizetmb:= (tmb/each_median)+0.1]
#   #each_tumor_info <- each_tumor_info[, normalizetmb:= tmb]
#   total_tumor_normalize <- rbindlist(list(total_tumor_normalize,each_tumor_info))
# }
# 
# total_cnv <- total_tumor_normalize

## Normalize end 

candidate_gene <- unique(total_tumor_gene_sum$gene_name)
# each tumor type seperate tmbl and tmbh
total_gene_sum <- NULL
for (each_gene in candidate_gene){
  #print(each_gene)
  tmb_l_group <- total_cnv[Subtype %in%tmb_l, ]
  tmb_l_group_total_samples <- unique(tmb_l_group$sample_names)
  tmb_l_gene_mut_samples <- unique(tmb_l_group[gene_name==each_gene,.(sample_names)][["sample_names"]])
  tmb_l_gene_no_mut_samples <- unique(setdiff(tmb_l_group_total_samples,tmb_l_gene_mut_samples))
  tmb_l_gene_mut_tmb <- unique(total_cnv[sample_names %in% tmb_l_gene_mut_samples,.(sample_names,tmb)])[["tmb"]]
  tmb_l_gene_no_mut_tmb <- unique(total_cnv[sample_names %in% tmb_l_gene_no_mut_samples,.(sample_names,tmb)])[["tmb"]]
  
  
  if (length(tmb_l_gene_mut_samples)!=0 && length(tmb_l_gene_no_mut_samples)!=0){
    tmb_l_test <- wilcox.test(tmb_l_gene_mut_tmb,tmb_l_gene_no_mut_tmb)
    tmb_l_pvalue <- tmb_l_test$p.value
    median_tmb_l_gene_mut_tmb <- median(tmb_l_gene_mut_tmb)
    median_tmb_l_gene_no_mut_tmb <- median(tmb_l_gene_no_mut_tmb)
    tmb_l_fold_change = median_tmb_l_gene_mut_tmb/median_tmb_l_gene_no_mut_tmb
  }
  # else{
  #   tmb_l_pvalue <- "No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
  #   median_tmb_l_gene_mut_tmb <- "No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
  #   median_tmb_l_gene_no_mut_tmb <-"No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
  # }
  
  tmb_h_group <- total_cnv[Subtype %in%tmb_h, ]
  tmb_h_group_total_samples <- unique(tmb_h_group$sample_names)
  tmb_h_gene_mut_samples <- unique(tmb_h_group[gene_name==each_gene,.(sample_names)][["sample_names"]])
  tmb_h_gene_no_mut_samples <- setdiff(tmb_h_group_total_samples,tmb_h_gene_mut_samples)
  tmb_h_gene_mut_tmb <- unique(total_cnv[sample_names %in% tmb_h_gene_mut_samples,.(sample_names,tmb)])[["tmb"]]
  tmb_h_gene_no_mut_tmb <- unique(total_cnv[sample_names %in% tmb_h_gene_no_mut_samples,.(sample_names,tmb)])[["tmb"]]
  
  if (length(tmb_h_gene_mut_samples)!=0 && length(tmb_h_gene_no_mut_samples)!=0){
    
    tmb_h_test <- wilcox.test(tmb_h_gene_mut_tmb,tmb_h_gene_no_mut_tmb)
    tmb_h_pvalue <- tmb_h_test$p.value
    median_tmb_h_gene_mut_tmb <- median(tmb_h_gene_mut_tmb)
    median_tmb_h_gene_no_mut_tmb <- median(tmb_h_gene_no_mut_tmb)
    tmb_h_fold_change = median_tmb_h_gene_mut_tmb/median_tmb_h_gene_no_mut_tmb
  }
  
  # else{
  #   tmb_h_pvalue <-"No tmb_h_gene_mut_samples or tmb_h_gene_no_mut_samples"
  #   median_tmb_h_gene_mut_tmb <- "No tmb_h_gene_mut_samples or tmb_h_gene_no_mut_samples"
  #   median_tmb_h_gene_no_mut_tmb <-"No tmb_h_gene_mut_samples or tmb_h_gene_no_mut_samples"
  # }
  
  # if (length(tmb_h_gene_mut_samples)!=0 && length(tmb_h_gene_no_mut_samples)!=0 
  #     && length(tmb_l_gene_mut_samples)!=0 && length(tmb_l_gene_no_mut_samples)!=0){
  each_gene_sum <- list(gene = each_gene,
                        median_tmb_l_gene_mut_tmb=median_tmb_l_gene_mut_tmb,
                        median_tmb_l_gene_no_mut_tmb=median_tmb_l_gene_no_mut_tmb,
                        median_tmb_h_gene_mut_tmb= median_tmb_h_gene_mut_tmb,
                        median_tmb_h_gene_no_mut_tmb=median_tmb_h_gene_no_mut_tmb,
                        tmb_l_number_samples=length(tmb_l_gene_mut_samples),
                        tmb_l_gene_no_mut_samples = length(tmb_l_gene_no_mut_samples),
                        tmb_h_gene_number_samples = length(tmb_h_gene_mut_samples),
                        tmb_h_gene_no_mut_samples = length(tmb_h_gene_no_mut_samples),
                        tmb_l_fold_change = tmb_l_fold_change, 
                        tmb_h_fold_change = tmb_h_fold_change,
                        tmb_l_pvalue = tmb_l_pvalue,
                        tmb_h_pvalue = tmb_h_pvalue)
  
  
  total_gene_sum <- rbindlist(list(total_gene_sum,each_gene_sum))
}
total_gene_sum <- setDT(total_gene_sum)
total_gene_sum <- total_gene_sum[order(tmb_l_pvalue)]
total_gene_sum$tmb_l_BH_pvalue <-  p.adjust(total_gene_sum$tmb_l_pvalue, method = "BH")
total_gene_sum <- total_gene_sum[order(tmb_h_pvalue)]
total_gene_sum$tmb_h_BH_pvalue <- p.adjust(total_gene_sum$tmb_h_pvalue,method = "BH")


fwrite(total_gene_sum, file = paste(output_dir,"03_02","amp_not_normalize_p_value_candidate_gene_associated_TMB_03_02.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")



#################### Combine SNV CNV indel  ################################
## 03_02 only examine pathway tp53 and cell cycle

base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

output_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/Gene_association_TMB/"

pathway <- fread(paste(base_dir,"all_pathway.txt",sep = seperator), na.strings = "")
target_path_way <- pathway[,c("TP53","Cell cycle"), with =F]
target_pathway_gene <- NULL
for (each_pathway in colnames(target_path_way)){
  each_path_gene <- pathway[, each_pathway, with = F][[each_pathway]]
  each_path_clean_gene <- each_path_gene[!is.na(each_path_gene)]
  target_pathway_gene <- c(target_pathway_gene,each_path_clean_gene)
}


#target_pathway_gene <- fread(paste(base_dir,"target_pathway_total_genes.txt",sep = seperator), na.strings = "")

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 
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

## amp_delete
amp_delete <- fread(paste(base_dir,"CNV_Taifang_total_amp_delete_no_pseudo_subtype_02_18.txt",sep = seperator),
                    header = T,na.strings = "")

amp_delete <- amp_delete[!sample_names %in% exclude]
amp_delete <- amp_delete[!grepl("ENSCAFG",amp_delete[,.(gene_name)]$gene_name,ignore.case = T)]
CNV <- unique(amp_delete[,c("sample_names","gene_name","mut_type","subtype"),with =F])
colnames(CNV)<- c("sample_names","gene_name","status","Subtype")
CNV <- CNV[gene_name %in% target_pathway_gene, ]

## exclude s1 high
s1_data <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
exclude <- c(exclude, s1_high_sample)
##
total_snv_cnv <- rbindlist(list(SNV,indel_file, CNV))
total_snv_cnv <- total_snv_cnv[!sample_names %in% exclude & Subtype !="UCL",]

all_tumor_cut <- 0
signle_tumor_cut <- 5

all_tumor_type <- unique(total_snv_cnv$Subtype)
# might only have one sample has that
total_tumor_gene_sum <- NULL
for( index in 1:length(all_tumor_type)){
  print(paste("processing the",index,"tumor, with total tumors",length(all_tumor_type),sep = " " ))
  each_tumor <- all_tumor_type[index]
  each_tumor_gene_info <- unique(total_snv_cnv[Subtype == each_tumor, .(gene_name)][["gene_name"]])
  each_tumor_gene_candidate <- NULL
  target_gene_sample_number <- NULL
  target_gene_total_sample_number <- NULL
  for (gene_index in 1:length(each_tumor_gene_info)){
    each_gene <- each_tumor_gene_info[gene_index]
    # each_gene <- "PIK3CA"
    total_sample <- unique(total_snv_cnv[gene_name==each_gene, .(sample_names,Subtype)])
    total_sample_sum <- as.data.table(table(total_sample$Subtype))
    #total_sample_number <- sum(total_sample_sum$N)
    each_tumor_sample_sum <- total_sample_sum[which(total_sample_sum$V1==each_tumor)]$N
    
    if (each_tumor_sample_sum >= signle_tumor_cut && total_sample_number>= all_tumor_cut){
      candidate_gene <- each_gene
      target_gene_sample_number <- c(target_gene_sample_number,each_tumor_sample_sum)
     
      each_tumor_gene_candidate <- c(each_tumor_gene_candidate,candidate_gene)
    }
  }
  each_tumor_sum <-  data.table(Subtype = each_tumor,
                                gene_name = each_tumor_gene_candidate,
                                target_gene_sample_number = target_gene_sample_number,
                                target_gene_total_sample_number = target_gene_total_sample_number) 
  
  total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_tumor_sum),fill = T)
}
total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)
fwrite(total_tumor_gene_sum, file = paste(output_dir,"03_02","include_amp_candidate_gene_associated_TMB_03_02.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

## Use candidate gene to compare tmb

## exclude s1 high
# s1_data <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
# s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
# exclude <- c(exclude, s1_high_sample)

## excldue NA and UCL
tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")

# total_tumor_gene_sum <- fread(paste(output_dir,"03_02","Not_include_amp_candidate_gene_associated_TMB_03_02.txt",sep = seperator))
# 
total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)

## append TMB info
TMB_info <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
colnames(TMB_info)
total_snv_cnv$tmb <- match_vector_table(total_snv_cnv$sample_names, "combine_snv_indel_tmb", TMB_info, string_value = F)
total_tumor_type <- unique(total_snv_cnv$Subtype)

candidate_gene <- unique(total_tumor_gene_sum$gene_name)
# each tumor type seperate tmbl and tmbh
total_gene_sum <- NULL
for (each_gene in candidate_gene){
  #print(each_gene)
  tmb_l_group <- total_snv_cnv[Subtype %in%tmb_l, ]
  tmb_l_group_total_samples <- unique(tmb_l_group$sample_names)
  tmb_l_gene_mut_samples <- unique(tmb_l_group[gene_name==each_gene,.(sample_names)][["sample_names"]])
  tmb_l_gene_no_mut_samples <- unique(setdiff(tmb_l_group_total_samples,tmb_l_gene_mut_samples))
  tmb_l_gene_mut_tmb <- unique(total_snv_cnv[sample_names %in% tmb_l_gene_mut_samples,.(sample_names,tmb)])[["tmb"]]
  tmb_l_gene_no_mut_tmb <- unique(total_snv_cnv[sample_names %in% tmb_l_gene_no_mut_samples,.(sample_names,tmb)])[["tmb"]]
  
  
  if (length(tmb_l_gene_mut_samples)!=0 && length(tmb_l_gene_no_mut_samples)!=0){
    tmb_l_test <- wilcox.test(tmb_l_gene_mut_tmb,tmb_l_gene_no_mut_tmb)
    tmb_l_pvalue <- tmb_l_test$p.value
    median_tmb_l_gene_mut_tmb <- median(tmb_l_gene_mut_tmb)
    median_tmb_l_gene_no_mut_tmb <- median(tmb_l_gene_no_mut_tmb)
    tmb_l_fold_change = median_tmb_l_gene_mut_tmb/median_tmb_l_gene_no_mut_tmb
  }
  # else{
  #   tmb_l_pvalue <- "No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
  #   median_tmb_l_gene_mut_tmb <- "No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
  #   median_tmb_l_gene_no_mut_tmb <-"No tmb_l_gene_mut_samples or tmb_l_gene_no_mut_samples"
  # }
  
  tmb_h_group <- total_snv_cnv[Subtype %in%tmb_h, ]
  tmb_h_group_total_samples <- unique(tmb_h_group$sample_names)
  tmb_h_gene_mut_samples <- unique(tmb_h_group[gene_name==each_gene,.(sample_names)][["sample_names"]])
  tmb_h_gene_no_mut_samples <- setdiff(tmb_h_group_total_samples,tmb_h_gene_mut_samples)
  tmb_h_gene_mut_tmb <- unique(total_snv_cnv[sample_names %in% tmb_h_gene_mut_samples,.(sample_names,tmb)])[["tmb"]]
  tmb_h_gene_no_mut_tmb <- unique(total_snv_cnv[sample_names %in% tmb_h_gene_no_mut_samples,.(sample_names,tmb)])[["tmb"]]
  
  if (length(tmb_h_gene_mut_samples)!=0 && length(tmb_h_gene_no_mut_samples)!=0){
    
    tmb_h_test <- wilcox.test(tmb_h_gene_mut_tmb,tmb_h_gene_no_mut_tmb)
    tmb_h_pvalue <- tmb_h_test$p.value
    median_tmb_h_gene_mut_tmb <- median(tmb_h_gene_mut_tmb)
    median_tmb_h_gene_no_mut_tmb <- median(tmb_h_gene_no_mut_tmb)
    tmb_h_fold_change = median_tmb_h_gene_mut_tmb/median_tmb_h_gene_no_mut_tmb
  }
  
  # else{
  #   tmb_h_pvalue <-"No tmb_h_gene_mut_samples or tmb_h_gene_no_mut_samples"
  #   median_tmb_h_gene_mut_tmb <- "No tmb_h_gene_mut_samples or tmb_h_gene_no_mut_samples"
  #   median_tmb_h_gene_no_mut_tmb <-"No tmb_h_gene_mut_samples or tmb_h_gene_no_mut_samples"
  # }
  
  # if (length(tmb_h_gene_mut_samples)!=0 && length(tmb_h_gene_no_mut_samples)!=0 
  #     && length(tmb_l_gene_mut_samples)!=0 && length(tmb_l_gene_no_mut_samples)!=0){
  each_gene_sum <- list(gene = each_gene,
                        median_tmb_l_gene_mut_tmb=median_tmb_l_gene_mut_tmb,
                        median_tmb_l_gene_no_mut_tmb=median_tmb_l_gene_no_mut_tmb,
                        median_tmb_h_gene_mut_tmb= median_tmb_h_gene_mut_tmb,
                        median_tmb_h_gene_no_mut_tmb=median_tmb_h_gene_no_mut_tmb,
                        tmb_l_number_samples=length(tmb_l_gene_mut_samples),
                        tmb_l_gene_no_mut_samples = length(tmb_l_gene_no_mut_samples),
                        tmb_h_gene_number_samples = length(tmb_h_gene_mut_samples),
                        tmb_h_gene_no_mut_samples = length(tmb_h_gene_no_mut_samples),
                        tmb_l_fold_change = tmb_l_fold_change, 
                        tmb_h_fold_change = tmb_h_fold_change,
                        tmb_l_pvalue = tmb_l_pvalue,
                        tmb_h_pvalue = tmb_h_pvalue)
  
  
  total_gene_sum <- rbindlist(list(total_gene_sum,each_gene_sum))
}
total_gene_sum <- setDT(total_gene_sum)
total_gene_sum <- total_gene_sum[order(tmb_l_pvalue)]
total_gene_sum$tmb_l_BH_pvalue <-  p.adjust(total_gene_sum$tmb_l_pvalue, method = "BH")
total_gene_sum <- total_gene_sum[order(tmb_h_pvalue)]
total_gene_sum$tmb_h_BH_pvalue <- p.adjust(total_gene_sum$tmb_h_pvalue,method = "BH")

a = total_gene_sum[total_gene_sum$gene=="MDM2",]


fwrite(total_gene_sum, file = paste(output_dir,"03_02","include_amp_not_normalize_p_value_candidate_gene_associated_TMB_03_02.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")




