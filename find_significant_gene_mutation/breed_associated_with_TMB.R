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

fwrite(total_breed_sum, file = paste(output_dir,"03_19","breeds_Not_include_amp_candidate_gene_associated_TMB_03_19_summary.txt",sep = seperator),
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
fwrite(total_breed_sum, file = paste(output_dir,"03_19","golden_withS1_associated_TMB_03_19_summary.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")


## Breed pathway  copy code from P53_cell cyle gene_associated_tmb.R:
## only need Breed,Pathway	Mutated_samples	P_value	Fold_change	BH_P_value
## for breed analysis, we only use snv +indel, no amp
## extract breed 01 matrix and compare tmb
base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

output_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/Gene_association_TMB/"

pathway <- fread(paste(base_dir,"all_pathway_03_19.txt",sep = seperator), na.strings = "")
path_col <- colnames(pathway)
target_path_way <- pathway[,c("p53","Cell cycle"), with =F]
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

total_sample <- unique(total_snv_cnv$sample_names)

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
total_sum$breed <- match_vector_table(total_sum$sample_names, "final_breed_label", whole_wes_clean_breed_table, string_value = F)

## normalize tmb across different tumor type
total_tumor_type <- unique(total_sum$Subtype)
total_tumor_normalize <- NULL

for (each_tumor in total_tumor_type){
  each_tumor_info <- total_sum[Subtype==each_tumor,]
  each_median <- median(unique(total_sum[Subtype==each_tumor,.(sample_names,tmb)])[['tmb']])
  #median(total_sum[Subtype==each_tumor, .(sample_names,tmb)][['tmb']])
  # each_sd <- sd(total_sum[Subtype==each_tumor, .(tmb)][['tmb']])
  each_tumor_info <- each_tumor_info[, normalizetmb:= (tmb/each_median)]
  #each_tumor_info <- each_tumor_info[, normalizetmb:= tmb]
  total_tumor_normalize <- rbindlist(list(total_tumor_normalize,each_tumor_info))
}


total_sum <- total_tumor_normalize

target_pathway_col <- c("p53","Cell cycle")
target_breeds <- c("Boxer","Cocker Spaniel","Golden Retriever","Greyhound","Maltese",
                   "Poodle","Rottweiler","Schnauzer","Shih Tzu",
                   "Yorkshire Terrier")
total_target_breed_type <- target_breeds
total_breed_pathway_sum <- NULL

for (breed_index in 1:length(total_target_breed_type)){
  
  each_breed <- total_target_breed_type[breed_index]
  each_breed_info <- total_sum[breed==each_breed,]
  each_breed_total_pathway_summary <- NULL
  for (each_pathway_index in 1:length(target_pathway_col)){
    each_breed_each_pathway_summary <- NULL
    each_pathway <- target_pathway_col[each_pathway_index]
    each_breed_pathway_total_sample_number <- sum(each_breed_info[, each_pathway, with = F])
    #print(c(each_breed,each_pathway))
    if (each_breed_pathway_total_sample_number >= 0 ){
      gene_mut_tmb <- each_breed_info[get(each_pathway) == 1, .(normalizetmb)][['normalizetmb']]
      gene_no_mut_tmb <- each_breed_info[get(each_pathway) == 0, .(normalizetmb)][['normalizetmb']]
      
      if (length(gene_mut_tmb)>0 & length(gene_no_mut_tmb)>0){
        tmb_test <- wilcox.test(gene_mut_tmb,gene_no_mut_tmb)
        p_value <- tmb_test$p.value
        fold_change <- median(gene_mut_tmb)/median(gene_no_mut_tmb)
        each_breed_each_pathway_summary <- list(Breed= each_breed,
                                                Pathway = each_pathway,
                                                Mutated_samples = each_breed_pathway_total_sample_number,
                                                P_value=p_value,
                                                Fold_change= fold_change)
        
      }
      else if(length(gene_mut_tmb)==0) {
        each_breed_each_pathway_summary <- list(Breed= each_breed,
                                                Pathway = each_pathway,
                                                Mutated_samples = 0,
                                                P_value="No_Pvalue",
                                                Fold_change= "No sample Mut")
      }
    }
    each_breed_total_pathway_summary <- rbindlist(list(each_breed_total_pathway_summary,each_breed_each_pathway_summary))  
    
  }
  each_breed_total_pathway_summary <- setDT(each_breed_total_pathway_summary)
  each_breed_total_pathway_summary <- each_breed_total_pathway_summary[order(P_value)]
  
  total_breed_pathway_sum <- rbindlist(list(total_breed_pathway_sum,each_breed_total_pathway_summary))
}

#each_breed_total_pathway_summary$BH_P_value <-  p.adjust(each_breed_total_pathway_summary$P_value, method = "BH")
## adjust p value##
final_out <- NULL
target_pathway_col
for (each_pathway in target_pathway_col){
  each_pathway_info <- total_breed_pathway_sum[Pathway == each_pathway]
  each_pathway_info <- each_pathway_info[order(P_value)]
  each_pathway_info$BH_P_value <-  p.adjust(each_pathway_info$P_value, method = "BH")
  final_out <- rbindlist(list(final_out,each_pathway_info))
}


fwrite(final_out, file = paste(output_dir,"03_19","all_breeds_cell_cycle_tp53_pathway.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

### data for TMB breed ##
breed_pathway_tmb <- total_sum[breed %in% target_breeds,.(sample_names,breed,p53,`Cell cycle`,tmb,Subtype)]
breed_pathway_tmb[breed_pathway_tmb==0] <- "Not altered"
breed_pathway_tmb[breed_pathway_tmb==1] <- "Altered"
fwrite(breed_pathway_tmb, file = paste(output_dir,"03_19","TMB_all_breeds_cell_cycle_tp53_pathway.txt",sep = seperator),
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
# fwrite(total_breed_sum, file = paste(output_dir,"03_19","Poodle_withS1_associated_TMB_03_19_summary.txt",sep = seperator),
#        col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")
# "No breed provided and unable to do the breed-prediction"

############################ Breed assoicated TMB end ############################
