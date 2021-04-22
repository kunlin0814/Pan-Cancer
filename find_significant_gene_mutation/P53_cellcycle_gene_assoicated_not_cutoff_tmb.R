# This script will look at the target gene association with TMB with no cut-off (unlike with cutoff script , cutoff=5) 
library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")
gene_assoicate_tmb <- function(each_tumor_info,each_tumor_uniq_gene,signle_tumor_cut,total_sample_number){
  each_tumor_total_gene_summary <- NULL
  for (gene_index in 1:length(each_tumor_uniq_gene)){
    each_tumor_each_gene_summary <- NULL
    each_gene <- each_tumor_uniq_gene[gene_index]
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
        each_tumor_each_gene_summary <- list( Gene = candidate_gene,
                                              Mutated_samples = each_tumor_gene_total_sample_number,
                                              P_value=p_value,
                                              Fold_change= fold_change,
                                              TP53_mutual_P_value=tp53_fisher_p,
                                              TP53_Incl_Excl=tp53_relation_type,
                                              TP53_mutual_significant = tp53_relation_sign)
        
        
      }
      else if (length(gene_mut_tmb)==0) {
        p_value <- "no smaple mutated"
        fold_change <- "no smaple mutated"
        each_tumor_each_gene_summary <- list(Gene = candidate_gene,
                                             Mutated_samples = each_tumor_gene_total_sample_number,
                                             P_value=p_value,
                                             Fold_change= fold_change,
                                             TP53_mutual_P_value="no smaple mutated",
                                             TP53_Incl_Excl="no smaple mutated",
                                             TP53_mutual_significant = "no smaple mutated")
      }
    }
    each_tumor_total_gene_summary <- rbindlist(list(each_tumor_total_gene_summary,each_tumor_each_gene_summary))  
    
  }
  each_tumor_total_gene_summary <- setDT(each_tumor_total_gene_summary)
  each_tumor_total_gene_summary <- each_tumor_total_gene_summary[order(P_value)]
  each_tumor_total_gene_summary$BH_P_value <-  p.adjust(each_tumor_total_gene_summary$P_value, method = "BH")
  return(each_tumor_total_gene_summary)
}


base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

output_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/Gene_association_TMB"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/Gene_association_TMB/"

pathway <- fread(paste(base_dir,"all_pathway_04_14.txt",sep = seperator), na.strings = "")

target_path_way <- pathway[,c("p53","Cell cycle"), with =F]
target_pathway_gene <- NULL
for (each_pathway in colnames(target_path_way)){
  each_path_gene <- pathway[, each_pathway, with = F][[each_pathway]]
  each_path_clean_gene <- each_path_gene[!is.na(each_path_gene)]
  target_pathway_gene <- c(target_pathway_gene,each_path_clean_gene)
}
target_pathway_gene <- unique(target_pathway_gene)
target_pathway_gene <- c("CDKN2A","ATM","CDK4")


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
#breed <- match_vector_table(total_mut$sample_names,"final_breed_label",whole_wes_clean_breed_table)

## amp_delete
amp_delete <- fread(paste(base_dir,"CNV_exclude_failQC_fixed_OM_total_amp_delete_no_pseudo_subtype_04_11.txt",sep = seperator),
                    header = T,na.strings = "")

amp_delete <- amp_delete[!sample_names %in% exclude]
amp_delete <- amp_delete[!grepl("ENSCAFG",amp_delete[,.(gene_name)]$gene_name,ignore.case = T)]
CNV <- unique(amp_delete[,c("sample_names","gene_name","mut_type","subtype"),with =F])
colnames(CNV)<- c("sample_names","gene_name","status","Subtype")

## combine SNV CNV indel
total_snv_cnv <- rbindlist(list(SNV,indel_file, CNV))
total_snv_cnv <- total_snv_cnv[!sample_names %in% exclude & Subtype !="UCL",]

## append tmb info
## append TMB 
TMB_info <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/all_pan-tumor_tmb_04_06.txt")
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/all_pan-tumor_tmb_04_06.txt")
colnames(TMB_info)
colnames(TMB_info)[1] <- "sample_names"
total_snv_cnv$tmb <- match_vector_table(total_snv_cnv$sample_names, "total_tmb", TMB_info, string_value = F)

## exclude s1 high and UCL
s1_data <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/S1_high_low.txt")
s1_high_sample <- s1_data[S1_Status=="S1 high"]$SampleName
exclude <- c(exclude, s1_high_sample)

## Normalize TMB with regards to each tumor median (in breed associated analysis)
total_tumor_type <- unique(total_snv_cnv$Subtype)
total_tumor_normalize <- NULL

for (each_tumor in total_tumor_type){
  each_tumor_info <- total_snv_cnv[Subtype==each_tumor,]
  each_median <- median(unique(total_snv_cnv[Subtype==each_tumor,.(sample_names,tmb)])[['tmb']])
  #median(total_mut[Subtype==each_tumor, .(sample_names,tmb)][['tmb']])
  # each_sd <- sd(total_mut[Subtype==each_tumor, .(tmb)][['tmb']])
  each_tumor_info <- each_tumor_info[, normalizetmb:= (tmb/each_median)]
  #each_tumor_info <- each_tumor_info[, normalizetmb:= tmb]
  total_tumor_normalize <- rbindlist(list(total_tumor_normalize,each_tumor_info))
}

total_snv_cnv <- total_tumor_normalize


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
## this one we use normalized tmb not regular tmb
for (gene_index in 1:length(pan_tumor_uniq_gene)){
  each_gene_summary <- list()
  each_gene <- pan_tumor_uniq_gene[gene_index]
  #each_gene <- "CDKN2B"
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

# fwrite(total_tumor_gene_sum, file = paste(output_dir,"04_14","include_amp_TP53_cell_cycle_gene_assication_TMB.txt",sep = seperator),
#        col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")
## cross tumor end

### examine each tumor type
total_snv_cnv[Subtype=="BCL" & gene_name=="CDK4",]

total_tumor_type <- unique(total_snv_cnv$Subtype)
compare_gene <- "TP53"
signle_tumor_cut <- 0
total_tumor_gene_sum <- NULL
total_sample_number <- length(unique(total_snv_cnv$sample_names))
for (tumor_index in 1:length(total_tumor_type)){
  each_tumor <- total_tumor_type[tumor_index]
  each_tumor_info <- total_snv_cnv[Subtype==each_tumor,]
  each_tumor_uniq_gene <- target_pathway_gene ## in this script we only exam two pathway gene
  each_tumor_total_gene_summary <- NULL
  each_tumor_sum <- gene_assoicate_tmb(each_tumor_info,each_tumor_uniq_gene,signle_tumor_cut,total_sample_number)
  each_tumor_sum[['Tumor_type']] <- each_tumor
  total_tumor_gene_sum <- rbindlist(list(total_tumor_gene_sum,each_tumor_sum))
}

fwrite(total_tumor_gene_sum, file = paste(output_dir,"04_14","test_not_cut_include_amp_TP53_cell_cycle_all_tumor_type_gene_assication_TMB.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

all_tumor_cut <- 0
signle_tumor_cut <- 0
tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")

tmb_l_unique_gene <- target_pathway_gene
tmb_l_group <- total_snv_cnv[Subtype %in% tmb_l,] 

### TMB-l candidate genes
total_tumor_gene_sum <- NULL
total_sample_number <- length(unique(tmb_l_group$sample_names))
tmbl_sum <- gene_assoicate_tmb(tmb_l_group,tmb_l_unique_gene,signle_tumor_cut,total_sample_number)

#total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)
fwrite(tmbl_sum, file = paste(output_dir,"04_14","test_not_cut_include_amp_tmb_l_candidate_gene_associated_TMB_04_22_summary.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")

### TMB-h 
tmb_h_group <- total_snv_cnv[Subtype %in% tmb_h,] 
total_tumor_gene_sum <- NULL
tmb_h_unique_gene <- target_pathway_gene

total_sample_number <- length(unique(tmb_h_group$sample_names))
tmbh_sum <- gene_assoicate_tmb(tmb_h_group,tmb_h_unique_gene,signle_tumor_cut,total_sample_number)

#total_tumor_gene_sum <- na.omit(total_tumor_gene_sum)
fwrite(tmbh_sum, file = paste(output_dir,"04_14","test_not_cut_include_amp_tmb_h_candidate_gene_associated_TMB_04_22_summary.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t", eol = "\n",na = "NA")


### test ###
tmb_l_group[gene_name =="CDKN2A"]
tmb_h_group[gene_name =="CDKN2A"]
total_snv_cnv[Subtype=="GLM" & gene_name =="ATM"]

