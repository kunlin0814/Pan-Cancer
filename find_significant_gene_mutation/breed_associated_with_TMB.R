library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")

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

## append TMB 
TMB_info <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
total_mut$tmb <- match_vector_table(total_mut$sample_names, "combine_snv_indel_tmb", TMB_info, string_value = F)

### append breeds information ###
breed <- match_vector_table(total_mut$sample_names,"final_breed_label",whole_wes_clean_breed_table,string_value = T)
total_mut$breeds <- breed

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

##  exclude UCL samples only but not exclude s1 high ##
exclude <- c(exclude)

withS1_total_mut <- total_mut[!sample_names %in% exclude & Subtype !="UCL",]

## exclude s1 high and UCL samples
s1_data <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
s1_exclude <- c(exclude, s1_high_sample)

withoutS1_total_mut <- total_mut[!sample_names %in% s1_exclude & Subtype !="UCL",]

############################ Breed assoicated TMB ############################

## select 5 sample of all genes
compare_gene <- "TP53"
signle_tumor_cut <- 5
breed_cut_off <- 3
pan_tumor_uniq_gene <- unique(withoutS1_total_mut$gene_name)
total_sample_number <- length(unique(withoutS1_total_mut$sample_names))
candidate_gene_list <- c()
for (gene_index in 1:length(pan_tumor_uniq_gene)){
  each_gene_summary <- list()
  each_gene <- pan_tumor_uniq_gene[gene_index]
  #print(each_gene)
  each_gene_total_sample <- unique(withoutS1_total_mut[gene_name==each_gene, .(sample_names)])
  each_gene_total_sample_number <- nrow(each_gene_total_sample)
  if (each_gene_total_sample_number >= signle_tumor_cut ){
    candidate_gene_list <- c(candidate_gene_list,each_gene)
  }
}
## each gene need 3 samples in each breed
target_breeds <- c("Boxer","Cocker Spaniel","Golden Retriever","Greyhound","Maltese",
                   "Poodle","Rottweiler","Schnauzer","Shih Tzu",
                   "Yorkshire Terrier")
candidate_gene_list <- unique(candidate_gene_list)
total_breed_sum <- NULL
for (each_candidate_breed in target_breeds){
  each_breed_sum <- NULL
  for (each_gene in candidate_gene_list){
    
    each_breed_each_gene <- NULL  
    each_breed_number <- nrow(unique(withoutS1_total_mut[gene_name==each_gene & breeds==each_candidate_breed,.(sample_names)]))
    
    if (each_breed_number >=breed_cut_off){
      
      candidate_gene <- each_gene
      target_breed_total_sample <- unique(withoutS1_total_mut[breeds==each_candidate_breed, ]$sample_names)
      target_breed_sample_with_mut <- unique(withoutS1_total_mut[gene_name==candidate_gene & breeds==each_candidate_breed,]$sample_names)
      target_breed_with_mut_tmb <- unique(withoutS1_total_mut[gene_name==candidate_gene & breeds==each_candidate_breed,.(sample_names,normalizetmb)])[['normalizetmb']]
      target_breed_sample_without_mut <- target_breed_total_sample[!target_breed_total_sample %in% target_breed_sample_with_mut]
      target_breed_sample_without_mut_tmb <- unique(withoutS1_total_mut[breeds==each_candidate_breed & sample_names %in% target_breed_sample_without_mut,.(sample_names,normalizetmb)])[['normalizetmb']]
      
      if (length(target_breed_with_mut_tmb) >0 & length(target_breed_sample_without_mut_tmb)>0){
        #print (c(each_candidate_breed,candidate_gene))
        tmb_test <- wilcox.test(target_breed_with_mut_tmb,target_breed_sample_without_mut_tmb)
        p_value <- tmb_test$p.value
        fold_change <- median(target_breed_with_mut_tmb)/median(target_breed_sample_without_mut_tmb)
        
        compare_gene_sample <- unique(withoutS1_total_mut[gene_name==compare_gene & breeds==each_candidate_breed,][['sample_names']])
        candidate_gene_sample<- unique(withoutS1_total_mut[gene_name==candidate_gene & breeds==each_candidate_breed][['sample_names']])
        alt_alt <- length(intersect(compare_gene_sample,candidate_gene_sample)) #TP53 and another gene sample
        alt_no_alt <- length(setdiff(compare_gene_sample,candidate_gene_sample)) # TP53 but not other gene mut sample
        no_alt_alt <- length(setdiff(candidate_gene_sample,compare_gene_sample))# other gene mut sample but not tp53 mut sample
        no_alt_no_alt <- total_sample_number-alt_alt-alt_no_alt-no_alt_alt
        tp53_fisher_test <- fisher.test(rbind(c(alt_alt, alt_no_alt), c(no_alt_alt, no_alt_no_alt))); 
        tp53_fisher_p <- tp53_fisher_test[["p.value"]];
        tp53_relation_type <- ifelse(tp53_fisher_test[["estimate"]] > 1, "Inclusive", "Exclusive");
        tp53_relation_sign <- ifelse(tp53_fisher_p <= 0.05, "Yes", "No");
        each_breed_each_gene <- list(Breeds = each_candidate_breed,
                                     Gene = candidate_gene,
                                     Mutated_samples = length(target_breed_sample_with_mut),
                                     P_value=p_value,
                                     Fold_change= fold_change,
                                     TP53_mutual_P_value=tp53_fisher_p,
                                     TP53_Incl_Excl=tp53_relation_type,
                                     TP53_mutual_significant = tp53_relation_sign)
        
      }
    }
    each_breed_sum <- rbindlist(list(each_breed_sum, each_breed_each_gene))
    
  }
  each_breed_sum <- setDT(each_breed_sum)
  each_breed_sum <- each_breed_sum[order(P_value)]
  each_breed_sum$BH_pvalue <-  p.adjust(each_breed_sum$P_value, method = "BH")
  total_breed_sum <- rbindlist(list(total_breed_sum,each_breed_sum))
}

fwrite(total_breed_sum, file = paste(output_dir,"03_02","breeds_Not_include_amp_candidate_gene_associated_TMB_03_02_summary.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

##### golden retriever with S1 ( not excldue S1) ####
target_breeds <- "Golden Retriever"
s1_golden <- withS1_total_mut[sample_names %in% s1_high_sample & breeds==target_breeds]
golden_candidate_gene_list <- unique(s1_golden$gene_name)
total_breed_sum <- NULL
for (each_candidate_breed in target_breeds){
  each_breed_sum <- NULL
  for (each_gene in golden_candidate_gene_list){
    
    each_breed_each_gene <- NULL  
    each_breed_number <- nrow(unique(s1_golden[gene_name==each_gene & breeds==each_candidate_breed,.(sample_names)]))
    
    if (each_breed_number >=breed_cut_off){
      
      candidate_gene <- each_gene
      target_breed_total_sample <- unique(s1_golden[breeds==each_candidate_breed, ]$sample_names)
      target_breed_sample_with_mut <- unique(s1_golden[gene_name==candidate_gene & breeds==each_candidate_breed,]$sample_names)
      target_breed_with_mut_tmb <- unique(s1_golden[gene_name==candidate_gene & breeds==each_candidate_breed,.(sample_names,normalizetmb)])[['normalizetmb']]
      target_breed_sample_without_mut <- target_breed_total_sample[!target_breed_total_sample %in% target_breed_sample_with_mut]
      target_breed_sample_without_mut_tmb <- unique(s1_golden[breeds==each_candidate_breed & sample_names %in% target_breed_sample_without_mut,.(sample_names,normalizetmb)])[['normalizetmb']]
      
      if (length(target_breed_with_mut_tmb) >0 & length(target_breed_sample_without_mut_tmb)>0){
        #print (c(each_candidate_breed,candidate_gene))
        tmb_test <- wilcox.test(target_breed_with_mut_tmb,target_breed_sample_without_mut_tmb)
        p_value <- tmb_test$p.value
        fold_change <- median(target_breed_with_mut_tmb)/median(target_breed_sample_without_mut_tmb)
        
        compare_gene_sample <- unique(s1_golden[gene_name==compare_gene & breeds==each_candidate_breed,][['sample_names']])
        candidate_gene_sample<- unique(s1_golden[gene_name==candidate_gene & breeds==each_candidate_breed][['sample_names']])
        alt_alt <- length(intersect(compare_gene_sample,candidate_gene_sample)) #TP53 and another gene sample
        alt_no_alt <- length(setdiff(compare_gene_sample,candidate_gene_sample)) # TP53 but not other gene mut sample
        no_alt_alt <- length(setdiff(candidate_gene_sample,compare_gene_sample))# other gene mut sample but not tp53 mut sample
        no_alt_no_alt <- total_sample_number-alt_alt-alt_no_alt-no_alt_alt
        tp53_fisher_test <- fisher.test(rbind(c(alt_alt, alt_no_alt), c(no_alt_alt, no_alt_no_alt))); 
        tp53_fisher_p <- tp53_fisher_test[["p.value"]];
        tp53_relation_type <- ifelse(tp53_fisher_test[["estimate"]] > 1, "Inclusive", "Exclusive");
        tp53_relation_sign <- ifelse(tp53_fisher_p <= 0.05, "Yes", "No");
        each_breed_each_gene <- list(Breeds = each_candidate_breed,
                                     Gene = candidate_gene,
                                     Mutated_samples = length(target_breed_sample_with_mut),
                                     P_value=p_value,
                                     Fold_change= fold_change,
                                     TP53_mutual_P_value=tp53_fisher_p,
                                     TP53_Incl_Excl=tp53_relation_type,
                                     TP53_mutual_significant = tp53_relation_sign)
        
      }
    }
    each_breed_sum <- rbindlist(list(each_breed_sum, each_breed_each_gene))
  }
  total_breed_sum <- rbindlist(list(total_breed_sum,each_breed_sum))
}
## Adjustment
total_breed_sum <- setDT(total_breed_sum)
total_breed_sum <- total_breed_sum[order(P_value)]
total_breed_sum$BH_pvalue <-  p.adjust(total_breed_sum$P_value, method = "BH")
fwrite(total_breed_sum, file = paste(output_dir,"03_02","golden_withS1_associated_TMB_03_02_summary.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

##### Poodle with S1 ( not excldue S1) no samples #######
# target_breeds <- "Poodle"
# s1_Poodle <- withS1_total_mut[sample_names %in% s1_high_sample & breeds==target_breeds]
# Poodle_candidate_gene_list <- unique(s1_Poodle$gene_name)
# total_breed_sum <- NULL
# for (each_candidate_breed in target_breeds){
#   each_breed_sum <- NULL
#   for (each_gene in golden_candidate_gene_list){
#     
#     each_breed_each_gene <- NULL  
#     each_breed_number <- nrow(unique(s1_golden[gene_name==each_gene & breeds==each_candidate_breed,.(sample_names)]))
#     
#     if (each_breed_number >=breed_cut_off){
#       
#       candidate_gene <- each_gene
#       target_breed_total_sample <- unique(s1_golden[breeds==each_candidate_breed, ]$sample_names)
#       target_breed_sample_with_mut <- unique(s1_golden[gene_name==candidate_gene & breeds==each_candidate_breed,]$sample_names)
#       target_breed_with_mut_tmb <- unique(s1_golden[gene_name==candidate_gene & breeds==each_candidate_breed,.(sample_names,normalizetmb)])[['normalizetmb']]
#       target_breed_sample_without_mut <- target_breed_total_sample[!target_breed_total_sample %in% target_breed_sample_with_mut]
#       target_breed_sample_without_mut_tmb <- unique(s1_golden[breeds==each_candidate_breed & sample_names %in% target_breed_sample_without_mut,.(sample_names,normalizetmb)])[['normalizetmb']]
#       
#       if (length(target_breed_with_mut_tmb) >0 & length(target_breed_sample_without_mut_tmb)>0){
#         #print (c(each_candidate_breed,candidate_gene))
#         tmb_test <- wilcox.test(target_breed_with_mut_tmb,target_breed_sample_without_mut_tmb)
#         p_value <- tmb_test$p.value
#         fold_change <- median(target_breed_with_mut_tmb)/median(target_breed_sample_without_mut_tmb)
#         
#         compare_gene_sample <- unique(s1_golden[gene_name==compare_gene & breeds==each_candidate_breed,][['sample_names']])
#         candidate_gene_sample<- unique(s1_golden[gene_name==candidate_gene & breeds==each_candidate_breed][['sample_names']])
#         alt_alt <- length(intersect(compare_gene_sample,candidate_gene_sample)) #TP53 and another gene sample
#         alt_no_alt <- length(setdiff(compare_gene_sample,candidate_gene_sample)) # TP53 but not other gene mut sample
#         no_alt_alt <- length(setdiff(candidate_gene_sample,compare_gene_sample))# other gene mut sample but not tp53 mut sample
#         no_alt_no_alt <- total_sample_number-alt_alt-alt_no_alt-no_alt_alt
#         tp53_fisher_test <- fisher.test(rbind(c(alt_alt, alt_no_alt), c(no_alt_alt, no_alt_no_alt))); 
#         tp53_fisher_p <- tp53_fisher_test[["p.value"]];
#         tp53_relation_type <- ifelse(tp53_fisher_test[["estimate"]] > 1, "Inclusive", "Exclusive");
#         tp53_relation_sign <- ifelse(tp53_fisher_p <= 0.05, "Yes", "No");
#         each_breed_each_gene <- list(Breeds = each_candidate_breed,
#                                      Gene = candidate_gene,
#                                      Mutated_samples = length(target_breed_sample_with_mut),
#                                      P_value=p_value,
#                                      Fold_change= fold_change,
#                                      TP53_mutual_P_value=tp53_fisher_p,
#                                      TP53_Incl_Excl=tp53_relation_type,
#                                      TP53_mutual_significant = tp53_relation_sign)
#         
#       }
#     }
#     each_breed_sum <- rbindlist(list(each_breed_sum, each_breed_each_gene))
#   }
#   total_breed_sum <- rbindlist(list(total_breed_sum,each_breed_sum))
# }
# ## Adjustment
# total_breed_sum <- setDT(total_breed_sum)
# total_breed_sum <- total_breed_sum[order(P_value)]
# total_breed_sum$BH_pvalue <-  p.adjust(total_breed_sum$P_value, method = "BH")
# fwrite(total_breed_sum, file = paste(output_dir,"03_02","Poodle_withS1_associated_TMB_03_02_summary.txt",sep = seperator),
#        col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")
# "No breed provided and unable to do the breed-prediction"

############################ Breed assoicated TMB end ############################
