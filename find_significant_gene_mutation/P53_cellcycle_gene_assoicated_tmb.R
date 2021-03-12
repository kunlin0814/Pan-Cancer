library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")

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
target_pathway_gene <- unique(target_pathway_gene)
whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 
exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

mutect_after_vaf <- fread(paste(base_dir,"NonSyn_Burair_filtering3_WithBreeds_Subtypes_QCpass_mutect_after_vaf_02_11.txt",
                                sep =seperator))

SNV <- unique(mutect_after_vaf[,c("sample_names","gene_name","status","Subtype"), with =F])

## indel
indel_file <- fread(paste(base_dir,"passQC_pan-tumor-total_indel_info_0214.txt",sep =seperator))
setcolorder(indel_file,c("sample_names","gene_name","emsembl_id","status"))
indel_file <- indel_file[,emsembl_id:=NULL]

Subtype <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table )
indel_file$Subtype <- Subtype


## amp_delete
amp_delete <- fread(paste(base_dir,"CNV_Taifang_total_amp_delete_no_pseudo_subtype_02_18.txt",sep = seperator),
                    header = T,na.strings = "")

amp_delete <- amp_delete[!sample_names %in% exclude]
amp_delete <- amp_delete[!grepl("ENSCAFG",amp_delete[,.(gene_name)]$gene_name,ignore.case = T)]
CNV <- unique(amp_delete[,c("sample_names","gene_name","mut_type","subtype"),with =F])
colnames(CNV)<- c("sample_names","gene_name","status","Subtype")


## exclude s1 high and UCL
s1_data <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
exclude <- c(exclude, s1_high_sample)

## combine SNV CNV indel
total_snv_cnv <- rbindlist(list(SNV,indel_file, CNV))
total_snv_cnv <- total_snv_cnv[!sample_names %in% exclude & Subtype !="UCL",]

## append tmb info

TMB_info <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
total_snv_cnv$tmb <- match_vector_table(total_snv_cnv$sample_names, "combine_snv_indel_tmb", TMB_info, string_value = F)

## Normalize TMB with regards to each tumor median (in breed associated analysis)
total_tumor_type <- unique(total_mut$Subtype)
total_tumor_normalize <- NULL

for (each_tumor in total_tumor_type){
  each_tumor_info <- total_mut[Subtype==each_tumor,]
  each_median <- median(unique(total_mut[Subtype==each_tumor,.(sample_names,tmb)])[['tmb']])
  #median(total_mut[Subtype==each_tumor, .(sample_names,tmb)][['tmb']])
  # each_sd <- sd(total_mut[Subtype==each_tumor, .(tmb)][['tmb']])
  each_tumor_info <- each_tumor_info[, normalizetmb:= (tmb/each_median)]
  #each_tumor_info <- each_tumor_info[, normalizetmb:= tmb]
  total_tumor_normalize <- rbindlist(list(total_tumor_normalize,each_tumor_info))
}

total_mut <- total_tumor_normalize



#################### Main code ########################
## Normalize end 
## need to normalize with median for each tumor in pan tumor analysis
## In this script , we only exam TP53 and cell cycle gene
## Pan_tumor all tumor types
compare_gene <- "TP53"
signle_tumor_cut <- 0
total_tumor_gene_sum <- NULL
pan_tumor_uniq_gene <- target_pathway_gene
total_sample_number <- length(unique(total_snv_cnv$sample_names))
for (gene_index in 1:length(pan_tumor_uniq_gene)){
  each_gene_summary <- list()
  each_gene <- pan_tumor_uniq_gene[gene_index]
  #each_gene <- "PIK3CA"
  each_gene_total_sample <- unique(total_snv_cnv[gene_name==each_gene, .(sample_names)])
  each_gene_total_sample_number <- nrow(each_gene_total_sample)
  if (each_gene_total_sample_number >= signle_tumor_cut ){
    candidate_gene <- each_gene
    gene_mut_tmb <- unique(total_snv_cnv[gene_name==candidate_gene, .(sample_names,normalizetmb)])[['normalizetmb']]
    gene_mut_sample <- unique(total_snv_cnv[gene_name==candidate_gene, .(sample_names)])[['sample_names']]
    
    gene_no_mut_tmb <-  unique(total_snv_cnv[!sample_names %in% gene_mut_sample, .(sample_names,normalizetmb)])[['normalizetmb']]
    if (length(gene_mut_tmb)>0 & length(gene_no_mut_tmb)>0){
    
    tmb_test <- wilcox.test(gene_mut_tmb,gene_no_mut_tmb)
    p_value <- tmb_test$p.value
    fold_change <- median(gene_mut_tmb)/median(gene_no_mut_tmb)
    
    compare_gene_sample <- unique(total_snv_cnv[gene_name==compare_gene,][['sample_names']])
    candidate_gene_sample<- unique(total_snv_cnv[gene_name==candidate_gene][['sample_names']])
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
    
  }
  total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_gene_summary))
}
total_tumor_gene_sum <- setDT(total_tumor_gene_sum)
total_tumor_gene_sum <- total_tumor_gene_sum[order(P_value)]
total_tumor_gene_sum$BH_P_value <-  p.adjust(total_tumor_gene_sum$P_value, method = "BH")

fwrite(total_tumor_gene_sum, file = paste(output_dir,"03_02","include_amp_TP53_cell_cycle_gene_assication_TMB.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")
## cross tumor end

### examine each tumor type
total_tumor_type <- unique(total_snv_cnv$Subtype)
compare_gene <- "TP53"
signle_tumor_cut <- 5
total_tumor_gene_sum <- NULL
total_sample_number <- length(unique(total_snv_cnv$sample_names))
for (tumor_index in 1:length(total_tumor_type)){
  each_tumor <- total_tumor_type[tumor_index]
  each_tumor_info <- total_snv_cnv[Subtype==each_tumor,]
  each_tumor_uniq_gene <- target_pathway_gene ## in this script we only exam two pathway gene
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
      if (length(gene_mut_tmb)>0 & length(gene_no_mut_tmb)>0){
      
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
    }
    each_tumor_total_gene_summary <- rbindlist(list(each_tumor_total_gene_summary,each_tumor_each_gene_summary))  
    
  }
  each_tumor_total_gene_summary <- setDT(each_tumor_total_gene_summary)
  each_tumor_total_gene_summary <- each_tumor_total_gene_summary[order(P_value)]
  each_tumor_total_gene_summary$BH_P_value <-  p.adjust(each_tumor_total_gene_summary$P_value, method = "BH")
  total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_tumor_total_gene_summary))
}

fwrite(total_tumor_gene_sum, file = paste(output_dir,"03_02","include_amp_TP53_cell_cycle_all_tumor_type_gene_assication_TMB.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

############################ TMB-l and TMB-h gene associated tmb ############################

all_tumor_cut <- 0
signle_tumor_cut <- 0
tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")

tmb_l_unique_gene <- target_pathway_gene
tmb_l_group <- total_snv_cnv[Subtype %in% tmb_l,] 
### TMB-l candidate genes
total_tumor_gene_sum <- NULL
total_sample_number <- length(unique(tmb_l_group$sample_names))
for (gene_index in 1:length(tmb_l_unique_gene)){
  each_gene_summary <- list()
  each_gene <- tmb_l_unique_gene[gene_index]
  #print(each_gene)
  each_gene_total_sample <- unique(tmb_l_group[gene_name==each_gene, .(sample_names)])
  each_gene_total_sample_number <- nrow(each_gene_total_sample)
  if (each_gene_total_sample_number >= signle_tumor_cut ){
    candidate_gene <- each_gene
    gene_mut_tmb <- unique(tmb_l_group[gene_name==candidate_gene, .(sample_names,tmb)])[['tmb']]
    gene_mut_sample <- unique(tmb_l_group[gene_name==candidate_gene, .(sample_names)])[['sample_names']]
    
    gene_no_mut_tmb <-  unique(tmb_l_group[!sample_names %in% gene_mut_sample, .(sample_names,tmb)])[['tmb']]
    if (length(gene_mut_tmb)>0 & length(gene_no_mut_tmb)>0){
    tmb_test <- wilcox.test(gene_mut_tmb,gene_no_mut_tmb)
    p_value <- tmb_test$p.value
    fold_change <- median(gene_mut_tmb)/median(gene_no_mut_tmb)
    
    compare_gene_sample <- unique(tmb_l_group[gene_name==compare_gene,][['sample_names']])
    candidate_gene_sample<- unique(tmb_l_group[gene_name==candidate_gene][['sample_names']])
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
  }
  total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_gene_summary))
}
total_tumor_gene_sum <- setDT(total_tumor_gene_sum)
total_tumor_gene_sum <- total_tumor_gene_sum[order(P_value)]
total_tumor_gene_sum$BH_P_value <-  p.adjust(total_tumor_gene_sum$P_value, method = "BH")

total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)
fwrite(total_tumor_gene_sum, file = paste(output_dir,"03_02","include_amp_tmb_l_candidate_gene_associated_TMB_03_02_summary.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

### TMB-h 
tmb_h_group <- total_snv_cnv[Subtype %in% tmb_h,] 
total_tumor_gene_sum <- NULL
tmb_h_unique_gene <- target_pathway_gene
total_sample_number <- length(unique(tmb_h_group$sample_names))
for (gene_index in 1:length(tmb_h_unique_gene)){
  each_gene_summary <- list()
  each_gene <- tmb_h_unique_gene[gene_index]
  #print(each_gene)
  each_gene_total_sample <- unique(tmb_h_group[gene_name==each_gene, .(sample_names)])
  each_gene_total_sample_number <- nrow(each_gene_total_sample)
  if (each_gene_total_sample_number >= signle_tumor_cut ){
    candidate_gene <- each_gene
    gene_mut_tmb <- unique(tmb_h_group[gene_name==candidate_gene, .(sample_names,tmb)])[['tmb']]
    gene_mut_sample <- unique(tmb_h_group[gene_name==candidate_gene, .(sample_names)])[['sample_names']]
    
    gene_no_mut_tmb <-  unique(tmb_h_group[!sample_names %in% gene_mut_sample, .(sample_names,tmb)])[['tmb']]
    if (length(gene_mut_tmb)>0 & length(gene_no_mut_tmb)>0){
    tmb_test <- wilcox.test(gene_mut_tmb,gene_no_mut_tmb)
    p_value <- tmb_test$p.value
    fold_change <- median(gene_mut_tmb)/median(gene_no_mut_tmb)
    
    compare_gene_sample <- unique(tmb_h_group[gene_name==compare_gene,][['sample_names']])
    candidate_gene_sample<- unique(tmb_h_group[gene_name==candidate_gene][['sample_names']])
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
  }
  total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_gene_summary))
}
total_tumor_gene_sum <- setDT(total_tumor_gene_sum)
total_tumor_gene_sum <- total_tumor_gene_sum[order(P_value)]
total_tumor_gene_sum$BH_P_value <-  p.adjust(total_tumor_gene_sum$P_value, method = "BH")

total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)
fwrite(total_tumor_gene_sum, file = paste(output_dir,"03_02","include_amp_tmb_h_candidate_gene_associated_TMB_03_02_summary.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

### Pathway whole TP53 and cell cycle
library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")

base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 
#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt")

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

## pathway data
pathway <- fread(paste(base_dir,"all_pathway.txt",sep = seperator), na.strings = "")
path_col <- colnames(pathway)

total_sample <- unique(total_snv_cnv$sample_names)

# create a total pathway list (0 and 1 in the list for each sample)
total_sum <- list()
for (i in path_col){
  total_sum[[i]] =numeric(length(total_sample))  
}
total_sum[["sample_names"]] =  character(length(total_sample))


for (index in 1:length(total_sample)){
  each_sample = total_sample[index]
  total_sum[["sample_names"]][index] = each_sample
  each_sample_gene <- total_snv_cnv[sample_names==each_sample]$gene_name
  
  for ( col_index in 1:length(path_col)){
    col_name <- path_col[col_index]
    each_path_way_all_gene <- pathway[, col_name, with = F][[col_name]]
    each_path_way_gene <- each_path_way_all_gene[!is.na(each_path_way_all_gene)]
    check_inside <- sum(each_sample_gene %in% each_path_way_gene)
    if (check_inside >0){
      total_sum[[col_name]][index] = 1
    }
    else{
      total_sum[[col_name]][index] = 0
    }
  }
}
# total_snv_cnv$TP53_pathway <- match_vector_table(total_snv_cnv$sample_names,"TP53",total_sum, string_value = F)
# total_snv_cnv$Cell_cycle_pathway <- match_vector_table(total_snv_cnv$sample_names,"Cell cycle",total_sum, string_value = F)

total_sum<- setDT(total_sum)
total_sum$Subtype <- match_vector_table(total_sum$sample_name,"DiseaseAcronym2", whole_wes_clean_breed_table)
#total_sum$Subtype <- subtype
total_sum$tmb <- match_vector_table(total_sum$sample_names, "combine_snv_indel_tmb", TMB_info, string_value = F)


### all tumor  ###
target_pathway_col <- c("TP53","Cell cycle")
total_tumor_type <- unique(total_sum$Subtype)
total_tumor_gene_sum <- NULL
total_sample_number <- length(unique(total_sum$sample_names))
for (tumor_index in 1:length(total_tumor_type)){
  each_tumor <- total_tumor_type[tumor_index]
  each_tumor_info <- total_sum[Subtype==each_tumor,]
  each_tumor_total_pathway_summary <- NULL
  for (each_pathway_index in 1:length(target_pathway_col)){
    each_tumor_each_pathway <- NULL
    each_pathway <- target_pathway_col[each_pathway_index]
    each_tumor_gene_total_sample_number <- sum(each_tumor_info[, each_pathway, with = F])
    if (each_tumor_gene_total_sample_number >= signle_tumor_cut ){
      gene_mut_tmb <- each_tumor_info[get(each_pathway) == 1, .(tmb)][['tmb']]
      gene_no_mut_tmb <- each_tumor_info[get(each_pathway) == 0, .(tmb)][['tmb']]

      if (length(gene_mut_tmb)>0 & length(gene_no_mut_tmb)>0){
        tmb_test <- wilcox.test(gene_mut_tmb,gene_no_mut_tmb)
        p_value <- tmb_test$p.value
        fold_change <- median(gene_mut_tmb)/median(gene_no_mut_tmb)
        each_tumor_each_gene_summary <- list(Tumor_type= each_tumor,
                                             Pathway = each_pathway,
                                             Mutated_samples = each_tumor_gene_total_sample_number,
                                             P_value=p_value,
                                             Fold_change= fold_change)
        
        
      }
    }
    each_tumor_total_pathway_summary <- rbindlist(list(each_tumor_total_pathway_summary,each_tumor_each_gene_summary))  
    
  }
  each_tumor_total_pathway_summary <- setDT(each_tumor_total_pathway_summary)
  each_tumor_total_pathway_summary <- each_tumor_total_pathway_summary[order(P_value)]
  each_tumor_total_pathway_summary$BH_P_value <-  p.adjust(each_tumor_total_pathway_summary$P_value, method = "BH")
  total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_tumor_total_pathway_summary))
}


fwrite(total_tumor_gene_sum, file = paste(output_dir,"03_02","include_amp_TP53_cell_cycle_all_tumor_type_pathway_association_TMB.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")



### tmb-l tmb-h  ###

tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")

tmb_l_pathway_group <- total_sum[Subtype %in% tmb_l]
tmb_h_pathway_group <- total_sum[Subtype %in% tmb_h]
tmb_l_pathway_sum <- NULL
tmb_h_pathway_sum <- NULL
for (each_pathway_index in 1:length(target_pathway_col)){
  
  each_pathway <- target_pathway_col[each_pathway_index]
  tmb_l_mut <- tmb_l_pathway_group[get(each_pathway) == 1, .(tmb)][['tmb']]
  tmb_l_not_mut <- tmb_l_pathway_group[get(each_pathway) == 0, .(tmb)][['tmb']]
  tmb_h_mut <- tmb_h_pathway_group[get(each_pathway) == 1, .(tmb)][['tmb']]
  tmb_h_not_mut <- tmb_h_pathway_group[get(each_pathway) == 0, .(tmb)][['tmb']]  
  
  if (length(tmb_l_mut)>0 & length(tmb_l_not_mut)>0){
    tmb_l_test <- wilcox.test(tmb_l_mut,tmb_l_not_mut)
    tmb_l_p_value <- tmb_l_test$p.value
    tmb_l_fold_change <- median(tmb_l_mut)/median(tmb_l_not_mut)
    tmb_l_summary <- list(Tumor_type= "tmb_l",
                                         Pathway = each_pathway,
                                         Mutated_samples = length(tmb_l_mut),
                                         P_value=tmb_l_p_value,
                                         Fold_change= tmb_l_fold_change)
  
    
    }
    if (length(tmb_h_mut)>0 & length(tmb_h_not_mut)>0){
      tmb_h_test <- wilcox.test(tmb_h_mut,tmb_h_not_mut)
      tmb_h_p_value <- tmb_h_test$p.value
      tmb_h_fold_change <- median(tmb_h_mut)/median(tmb_h_not_mut)
      tmb_h_summary <- list(Tumor_type= "tmb_h",
                                           Pathway = each_pathway,
                                           Mutated_samples = length(tmb_h_mut),
                                           P_value=tmb_h_p_value,
                                           Fold_change= tmb_h_fold_change)
    }
  
  tmb_l_pathway_sum <- rbindlist(list(tmb_l_pathway_sum, tmb_l_summary))
  tmb_h_pathway_sum <- rbindlist(list(tmb_h_pathway_sum, tmb_h_summary))
  
}
tmb_l_pathway_sum <- setDT(tmb_l_pathway_sum)
tmb_l_pathway_sum <- tmb_l_pathway_sum[order(P_value)]
tmb_l_pathway_sum$BH_P_value <-  p.adjust(tmb_l_pathway_sum$P_value, method = "BH")

fwrite(tmb_l_pathway_sum, file = paste(output_dir,"03_02","include_amp_TP53_cell_cycle_tmb_l_type_pathway_association_TMB.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")
tmb_h_pathway_sum <- setDT(tmb_h_pathway_sum)
tmb_h_pathway_sum <- tmb_h_pathway_sum[order(P_value)]
tmb_h_pathway_sum$BH_P_value <-  p.adjust(tmb_h_pathway_sum$P_value, method = "BH")
fwrite(tmb_h_pathway_sum, file = paste(output_dir,"03_02","include_amp_TP53_cell_cycle_tmb_h_type_pathway_association_TMB.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")
