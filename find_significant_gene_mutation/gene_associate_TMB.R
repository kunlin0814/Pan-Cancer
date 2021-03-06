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

whole_wes_clean_breed_table <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt")
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt") 
exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

mutect_after_vaf <- fread(paste(base_dir,"Final_Total_withGene_Burair_Filtering3_VAF_Mutect_orientBiasModified_04_02.txt.gz",
                                sep =seperator))

a = unique(whole_wes_clean_breed_table[The_reason_to_exclude=="Pass QC" & final_breed_label=="Golden Retriever",.(Case_ID,final_breed_label)])
## append columns
mutect_after_vaf <- mutect_after_vaf[status!= "synonymous",]
mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude]
mutect_after_vaf <- mutect_after_vaf[tumor_type!="UCL"]
mutect_after_vaf$Subtype <- match_vector_table(mutect_after_vaf$sample_names,column = "DiseaseAcronym2", table =whole_wes_clean_breed_table,string_value = T )
mutect_after_vaf$finalbreed <- match_vector_table(mutect_after_vaf$sample_names,column="final_breed_label",table=whole_wes_clean_breed_table)

## indel
indel_file <- fread(paste(base_dir,"total_CDS_indel_info_withGene_04_08.txt",sep =seperator))
colnames(indel_file) <- c('chrom','pos','ref','alt','gene_name','ensembl_id','status','sample_names')
indel_file <- indel_file[!sample_names %in% exclude,]
#indel_file <- indel_file[gene_name!="-" & status!="nonframeshift " & ! sample_names %in% exclude,]
indel_file$Subtype <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table)
indel_file$finalbreed <- match_vector_table(indel_file$sample_names,"final_breed_label",whole_wes_clean_breed_table)
indel_file <- indel_file[gene_name!="-" & status=="frameshift" & ! sample_names %in% exclude,]
indel_file <- indel_file[,c("sample_names","gene_name","status","Subtype"),with=F]
#setcolorder(indel_file,c("sample_names","gene_name","emsembl_id","status"))
#indel_file <- indel_file[,emsembl_id:=NULL]

SNV <- unique(mutect_after_vaf[,c("sample_names","gene_name","status","Subtype"), with =F])
total_mut <- rbindlist(list(SNV,indel_file))

total_mut <- total_mut[!sample_names %in% exclude & Subtype !="UCL",]
breed <- match_vector_table(total_mut$sample_names,"final_breed_label",whole_wes_clean_breed_table)

total_mut$breed <- breed

a = unique(total_mut[Subtype=="MT",.(sample_names,breed)])
a <- setDT(a)
b = a[,.N,keyby= .(breed)]

unique(total_mut[gene_name=="TP53" & Subtype=='OSA',.(sample_names)])
length(unique(total_mut[Subtype=="OSA"]$sample_names))

## exclude s1 high and UCL samples
s1_data <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
a = match_vector_table(s1_high_sample, "DiseaseAcronym2",whole_wes_clean_breed_table)
b = data.frame(subtype = a, sample = s1_high_sample)
exclude <- c(exclude, s1_high_sample)
total_mut <- total_mut[!sample_names %in% exclude & Subtype !="UCL",]
## append TMB 
TMB_info <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/all_pan-tumor_tmb_04_06.txt")
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/all_pan-tumor_tmb_04_06.txt")
colnames(TMB_info)
colnames(TMB_info)[1] <- "sample_names"
total_mut$tmb <- match_vector_table(total_mut$sample_names, "total_tmb", TMB_info, string_value = F)

## Normalize TMB with regards to each tumor median (in pan-tumor analysis)
## now decide not normalize 2/25
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

#check result#

### test ###
# tar_gene_name <- "TP53"
# tumor_type <- "OSA"
# tar_info <- total_mut[Subtype==tumor_type]
# a = tar_info[sample_names=="GarLar",]
# mut_sample <- unique(tar_info [gene_name == tar_gene_name,.(sample_names)])$sample_names
# mut_sample_tmb <- unique(tar_info[sample_names %in% mut_sample,.(sample_names,tmb)])$tmb
# 
# not_mut_samples <- unique(tar_info[!sample_names %in% mut_sample,.(sample_names)])
# not_mut_sample_tmb <- unique(tar_info[!sample_names %in% mut_sample,.(sample_names,tmb)])$tmb
# median(mut_sample_tmb)/median(not_mut_sample_tmb)
# 
# wilcox.test(mut_sample_tmb,not_mut_sample_tmb)
# 
# write.table(mut_sample,file = 'C:/Users/abc73/Desktop/OSA_mut.txt',
#        col.names = F, row.names = F, quote = F, sep = "\n")
# 
# tp53 <- unique(total_mut[gene_name =="TP53",.(sample_names,)])
# tp53_sample <- tp53$sample_names
# no_tp53 <- setdiff(unique(total_mut$sample_names), tp53_sample)
# tp53_tmb <- tp53$normalizetmb
# no_tp53_tmb <- unique(total_mut[sample_names %in%no_tp53, .(sample_names,normalizetmb)])
# 
# wilcox.test(tp53_tmb,no_tp53_tmb$normalizetmb)
# median(tp53_tmb)/median(no_tp53_tmb$normalizetmb)

##

## Normalize end 


## identify candidate genes (5 samples across tumor type)

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
    #candidate_gene <- 'TP53'
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

fwrite(total_tumor_gene_sum, file = paste(output_dir,"04_14","pan_tumor_gene_assication_TMB.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")
## cross tumor end
nrow(unique(total_mut[gene_name=='PIK3CA' & Subtype=="HSA",.(sample_names)]))
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
 
fwrite(total_tumor_gene_sum, file = paste(output_dir,"04_14","all_tumor_type_gene_assication_TMB.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

############################ TMB-l and TMB-h gene associated tmb ############################

## Gene assocaited TMB only exam somatic mutation
base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

output_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/Gene_association_TMB"
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/Gene_association_TMB/"

whole_wes_clean_breed_table <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt")
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt") 
exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

mutect_after_vaf <- fread(paste(base_dir,"Final_Total_withGene_Burair_Filtering3_VAF_Mutect_orientBiasModified_04_02.txt.gz",
                                sep =seperator))
## append columns
mutect_after_vaf <- mutect_after_vaf[status!= "synonymous",]
mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude]
mutect_after_vaf <- mutect_after_vaf[tumor_type!="UCL"]
mutect_after_vaf$Subtype <- match_vector_table(mutect_after_vaf$sample_names,column = "DiseaseAcronym2", table =whole_wes_clean_breed_table,string_value = T )
mutect_after_vaf$finalbreed <- match_vector_table(mutect_after_vaf$sample_names,column="final_breed_label",table=whole_wes_clean_breed_table)
a = total_mut[,.N,keyby=.(Subtype,gene_name)][order(Subtype,-N)]
## indel
indel_file <- fread(paste(base_dir,"total_CDS_indel_info_withGene_04_08.txt",sep =seperator))
colnames(indel_file) <- c('chrom','pos','ref','alt','gene_name','ensembl_id','status','sample_names')
indel_file <- indel_file[!sample_names %in% exclude,]
#indel_file <- indel_file[gene_name!="-" & status!="nonframeshift " & ! sample_names %in% exclude,]
indel_file$Subtype <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table)
indel_file$finalbreed <- match_vector_table(indel_file$sample_names,"final_breed_label",whole_wes_clean_breed_table)
indel_file <- indel_file[gene_name!="-" & status!="nonframeshift" & ! sample_names %in% exclude,]
indel_file <- indel_file[,c("sample_names","gene_name","status","Subtype"),with=F]
#setcolorder(indel_file,c("sample_names","gene_name","emsembl_id","status"))
#indel_file <- indel_file[,emsembl_id:=NULL]

SNV <- unique(mutect_after_vaf[,c("sample_names","gene_name","status","Subtype"), with =F])
total_mut <- rbindlist(list(SNV,indel_file))

total_mut <- total_mut[!sample_names %in% exclude & Subtype !="UCL",]
breed <- match_vector_table(total_mut$sample_names,"final_breed_label",whole_wes_clean_breed_table)
total_mut$breed <- breed
## exclude s1 high and UCL samples
s1_data <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
exclude <- c(exclude, s1_high_sample)

total_mut <- total_mut[!sample_names %in% exclude & Subtype !="UCL",]
## append TMB 
TMB_info <- fread('G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/all_pan-tumor_tmb_04_06.txt')
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/all_pan-tumor_tmb_04_06.txt")
colnames(TMB_info)[1] <- "sample_names"
#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/all_pan-tumor_tmb_04_06.txt")
total_mut$tmb <- match_vector_table(total_mut$sample_names, "total_tmb", TMB_info, string_value = F)

all_tumor_cut <- 0
signle_tumor_cut <- 5
tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")

tmb_l_unique_gene <- unique(total_mut[Subtype %in% tmb_l,.(gene_name)][['gene_name']])
tmb_l_group <- total_mut[Subtype %in% tmb_l,] 
tmb_l_group[gene_name=='KAT6B',]
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
  total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_gene_summary))
}
total_tumor_gene_sum <- setDT(total_tumor_gene_sum)
total_tumor_gene_sum <- total_tumor_gene_sum[order(P_value)]
total_tumor_gene_sum$BH_P_value <-  p.adjust(total_tumor_gene_sum$P_value, method = "BH")

total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)
fwrite(total_tumor_gene_sum, file = paste(output_dir,"04_14","tmb_l_Not_include_amp_candidate_gene_associated_TMB_03_19_summary.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

### TMB-h 
tmb_h_group <- total_mut[Subtype %in% tmb_h,] 
total_tumor_gene_sum <- NULL
tmb_h_unique_gene <- unique(tmb_h_group$gene_name)
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
  total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_gene_summary))
}
total_tumor_gene_sum <- setDT(total_tumor_gene_sum)
total_tumor_gene_sum <- total_tumor_gene_sum[order(P_value)]
total_tumor_gene_sum$BH_P_value <-  p.adjust(total_tumor_gene_sum$P_value, method = "BH")

total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)
fwrite(total_tumor_gene_sum, file = paste(output_dir,"04_14","tmb_h_Not_include_amp_candidate_gene_associated_TMB_03_19_summary.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

## check results 
test_gene = "TP53"
tumor = "MT"
mut_sample <- unique(total_mut[gene_name==test_gene & Subtype==tumor]$sample_names)
not_mut_sample <- unique(total_mut[!sample_names %in% mut_sample & Subtype==tumor]$sample_names)
mut_tmb <- unique(total_mut[sample_names %in% mut_sample,.(sample_names,tmb)])
not_mut_tmb <- unique(total_mut[sample_names %in% not_mut_sample,.(sample_names,tmb)])
median(mut_tmb$tmb)/median(not_mut_tmb$tmb)
wilcox.test(mut_tmb$tmb,not_mut_tmb$tmb )

## P53 and cell cycle write in a different script 03/02

### P53 and cell cycle gene 

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

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt") 
exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

mutect_after_vaf <- fread(paste(base_dir,"Final_Total_withGene_Burair_Filtering3_VAF_Mutect_orientBiasModified_04_02.txt.gz",
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

## exclude s1 high and UCL
s1_data <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
exclude <- c(exclude, s1_high_sample)
##
total_snv_cnv <- rbindlist(list(SNV,indel_file, CNV))
total_snv_cnv <- total_snv_cnv[!sample_names %in% exclude & Subtype !="UCL" & gene_name %in% target_pathway_gene,]

## append tmb info

TMB_info <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/all_pan-tumor_tmb_04_06.txt")
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/all_pan-tumor_tmb_04_06.txt")
total_snv_cnv$tmb <- match_vector_table(total_snv_cnv$sample_names, "total_tmb", TMB_info, string_value = F)

## Normalize TMB with regards to each tumor median (in breed associated analysis)
all_tumor_type <- unique(total_snv_cnv$Subtype)
total_tumor_normalize <- NULL

for (each_tumor in all_tumor_type){
  each_tumor_info <- total_snv_cnv[Subtype==each_tumor,]
  each_median <- median(total_snv_cnv[Subtype==each_tumor, .(tmb)][['tmb']])
  # each_sd <- sd(total_snv_cnv[Subtype==each_tumor, .(tmb)][['tmb']])
  each_tumor_info <- each_tumor_info[, normalizetmb:= (tmb/each_median)+0.1]
  #each_tumor_info <- each_tumor_info[, normalizetmb:= tmb]
  total_tumor_normalize <- rbindlist(list(total_tumor_normalize,each_tumor_info))
}

total_snv_cnv <- total_tumor_normalize


all_tumor_cut <- 0
signle_tumor_cut <- 5

# all_tumor_cut <- 10
# signle_tumor_cut <- 5

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
fwrite(total_tumor_gene_sum, file = paste(output_dir,"05_01","include_amp_candidate_gene_associated_TMB_05_01.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

## Use candidate gene to compare tmb

## exclude s1 high
# s1_data <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
# s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
# exclude <- c(exclude, s1_high_sample)

## excldue NA and UCL
tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")

# total_tumor_gene_sum <- fread(paste(output_dir,"02_25","Not_include_amp_candidate_gene_associated_TMB_02_25.txt",sep = seperator))
# 
total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)

## append TMB info
TMB_info <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/all_pan-tumor_tmb_04_06.txt")
colnames(TMB_info)
total_cnv$tmb <- match_vector_table(total_cnv$sample_names, "total_tmb", TMB_info, string_value = F)
total_tumor_type <- unique(total_cnv$Subtype)

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


fwrite(total_gene_sum, file = paste(output_dir,"05_01","include_amp_not_normalize_p_value_candidate_gene_associated_TMB_05_01.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")


