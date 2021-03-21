# the script is for table and figure3 in pan-cancer
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
pathway <- fread(paste(base_dir,"all_pathway_03_19.txt",sep = seperator), na.strings = "")
path_col <- colnames(pathway)

## mutation Data

## SNV
mutect_after_vaf <- fread(paste(base_dir,"NonSyn_Burair_filtering3_WithBreeds_Subtypes_QCpass_mutect_after_vaf_02_11.txt",
                                sep =seperator))
## indel
  
indel_file <- fread(paste(base_dir,"passQC_pan-tumor-total_indel_info_0214.txt",sep =seperator))
indel_file <- indel_file[gene_name!="-" & status=="frameshift" & ! sample_names %in% exclude,]
setcolorder(indel_file,c("sample_names","gene_name","emsembl_id","status"))
indel_file <- indel_file[,emsembl_id:=NULL]
# fwrite(indel_file, file= paste(base_dir,"passQC_pan-tumor-total_indel_info_0214.txt",sep =seperator),
#        col.names = T, row.names = F, sep = "\t", quote = F, eol="\n")

Subtype <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table )
indel_file$Subtype <- Subtype

## CNV data
amp_delete <- fread(paste(base_dir,"CNV_Taifang_total_amp_delete_no_pseudo_subtype.txt",sep = seperator),header = T)
SNV <- unique(mutect_after_vaf[,c("sample_names","gene_name","status","Subtype"), with =F])
CNV <- unique(amp_delete[,c("sample_names","gene_name","mut_type","subtype"),with =F])
colnames(CNV)<- c("sample_names","gene_name","status","Subtype")


## combine snv, cnv, amp
total_mut <- rbindlist(list(SNV,indel_file,CNV))

total_sample <- unique(total_mut$sample_names)

total_mut <- total_mut[!sample_names %in% exclude,]

a = unique(total_mut[Subtype=="GLM",.(sample_names)])
b = unique(total_mut[Subtype=="HSA",.(sample_names)])

# create a total pathway list (0 and 1 in the list for each sample)
total_sum <- list()
for (i in path_col){
  total_sum[[i]] =numeric(length(total_sample))  
}
total_sum[["sample_names"]] =  character(length(total_sample))


#each_sample_gene <- total_mut[sample_names=="004"]$gene_name
for (index in 1:length(total_sample)){
  each_sample = total_sample[index]
  total_sum[["sample_names"]][index] = each_sample
  each_sample_gene <- total_mut[sample_names==each_sample]$gene_name
  
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

total_sum<- setDT(total_sum)
subtype <- match_vector_table(total_sum$sample_name,"DiseaseAcronym2", whole_wes_clean_breed_table)
total_sum$Subtype <- subtype

## exclude UCL tumor
no_UCL <- total_sum[Subtype!="UCL",]

## append TMB info
# check PIK3CA mut between tmb-h and tmb-l

TMB_info <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate/With_Breeds_exclude_failQC_TMB_Burair_filtering3_02_11.txt")
colnames(TMB_info)
no_UCL$tmb <- match_vector_table(no_UCL$sample_names, "combine_snv_indel_tmb", TMB_info, string_value = F)

no_UCL[no_UCL==1] <- "Alter"
no_UCL[no_UCL==0] <- "No Alter"
tmb_l <- c("MT", "GLM","BCL")
tmb_h <- c("OM","OSA","HSA","TCL")

insert = character(length(no_UCL$Subtype))
tmb_l_index <- which(no_UCL$Subtype %in% tmb_l)
tmb_h_index <- which(no_UCL$Subtype %in% tmb_h)
insert[tmb_l_index] <- "tmb_l"
insert[tmb_h_index] <- "tmb_h"
no_UCL$class <- insert


fwrite(no_UCL, file=paste(base_dir,"02_25","figureS6B_TMB_pathway_03_21.txt",sep=seperator),
       col.names = T, row.names = F, quote = F, sep="\t",eol="\n",na = "NA")

#tmb_l_group
# tmbl
tmb_l_group <- no_UCL[Subtype %in% tmb_l]

#withmut_l <- tmb_l_group[PIK3CA_Gene ==1,]$tmb
#withoutmut_l <- tmb_l_group[PIK3CA_Gene ==0,]$tmb
#wilcox.test(withmut_l,withoutmut_l)
#median(withmut_l)
#median(withoutmut_l)
# tmbh
tmb_h_group <- no_UCL[Subtype %in% tmb_h]
#withmut_h <- tmb_h_group[PIK3CA_Gene ==1,]$tmb
#withoutmut_h <- tmb_h_group[PIK3CA_Gene ==0,]$tmb
#wilcox.test(withmut_h,withoutmut_h)
#median(withmut_h)
#median(withoutmut_h)
# 
# fwrite(no_UCL,file = paste(base_dir,"02_25","no_UCL_pathway_summary_02_25.txt",sep = seperator),
#        col.names = T, row.names = F, quote = F, sep = "\t", na="NA")

total_tumor_p_value_sum <- NULL
total_tumor_type <- unique(no_UCL$Subtype)
for (index in 1:length(total_tumor_type)){
  
  each_tumor <- total_tumor_type[index]
  total_target_sample_number <- length(unique(no_UCL[Subtype==each_tumor]$sample_name))
  total_other_sample_number <- length(unique(no_UCL[Subtype!=each_tumor]$sample_name))
  each_tumor_target <- numeric(length(path_col))
  each_tumor_target_without <- numeric(length(path_col))
  each_tumor_other_with <- numeric(length(path_col))
  each_tumor_other_without <- numeric(length(path_col))
  each_tumor_enrich_scores <- numeric(length(path_col))
  p_value_arr <- numeric(length(path_col))
  odds_ratio_arr <- numeric(length(path_col))
  #each_tumor_enrich_p_value <- numeric(length(path_col))
  #each_tumor_deplete_p_value <- numeric(length(path_col))
  each_tumor_col <- character(length(path_col))
  each_tumor_type <- character(length(path_col))
  for (col_index in 1:length(path_col)){
    each_col <- path_col[col_index] 
    each_tumor_col[col_index] <- path_col[col_index]
    # consdier all 0 or all 1
    # use total and total- with = without
    # sum = each_tumor_target
    #total_sample <- sum(data.table(table(no_UCL[Subtype==each_tumor,each_col, with =F]))$N)
    each_tumor_type[col_index] <- each_tumor
    each_tumor_target[col_index] <- sum(no_UCL[Subtype==each_tumor,each_col, with =F])
    each_tumor_target_without[col_index] <- total_target_sample_number-each_tumor_target[col_index]
    each_tumor_other_with[col_index] <- sum(no_UCL[Subtype!=each_tumor,each_col, with =F])
    each_tumor_other_without[col_index] <- total_other_sample_number-each_tumor_other_with[col_index]
    
    testor <- rbind(c(each_tumor_target[col_index],each_tumor_target_without[col_index]),
                    c(each_tumor_other_with[col_index],each_tumor_other_without[col_index]))
    
    fisher_test <- fisher.test(testor)
    p_value <- fisher_test[["p.value"]];
    p_value_arr[col_index] <- p_value
    odds_ratio <- fisher_test[["estimate"]];
    odds_ratio_arr[col_index] <- odds_ratio
    enr_score <- min(5, (log10(p_value) * -1));
    if(odds_ratio < 1) {
      enr_score <- 0-enr_score;
    }
    # enrich_p_value <- fisher.test(testor, alternative = "greater")$p.value
    # deplete_p_value <- fisher.test(testor, alternative = "less")$p.value
    each_tumor_enrich_scores[col_index] <- enr_score
    # each_tumor_enrich_p_value[col_index] <- enrich_p_value
    # each_tumor_deplete_p_value[col_index] <- deplete_p_value
    #print(each_col)
  }
  each_tumor_sum <- data.table(tumor_type = each_tumor_type,
                               gene_pathway = each_tumor_col,
                               each_gene_with = each_tumor_target,
                               each_gene_without = each_tumor_target_without,
                               each_gene_other_with = each_tumor_other_with,
                               each_gene_other_without = each_tumor_other_without,
                               fisher_pvalue = p_value_arr,
                               odds_ratio = odds_ratio_arr,
                               each_tumor_enrich_scores= each_tumor_enrich_scores)
                               # each_gene_enrich_pvalue = each_tumor_enrich_p_value,
                               # each_gene_deplete_pvalue= each_tumor_deplete_p_value
  
  # 
  # each_tumor_sum <- each_tumor_sum[order(each_gene_enrich_pvalue)]
  # each_tumor_sum$Enrich_BH_pvalue = p.adjust(each_tumor_sum$each_gene_enrich_pvalue, method = "BH")
  # each_tumor_sum <- each_tumor_sum[order(each_gene_deplete_pvalue)]
  # each_tumor_sum$Deplete_BH_pvalue = p.adjust(each_tumor_sum$each_gene_deplete_pvalue, method = "BH")

  total_tumor_p_value_sum <- rbind(total_tumor_p_value_sum, each_tumor_sum)
  
}
#total_tumor_p_value_sum <- total_tumor_p_value_sum[gene_pathway %in% target_pathway,]

fwrite(total_tumor_p_value_sum, file=paste(base_dir,"02_25","With_pValue_pan-tumor_pathway_02_25.txt",sep=seperator),
       col.names = T, row.names = F, quote = F, sep="\t",eol="\n",na = "NA")

#### reshape and heatmap ###
# target_col <- path_col[-c(2,3,7,8)]
# target <- total_tumor_p_value_sum[gene_pathway %in% target_col,.(tumor_type,each_tumor_enrich_scores,gene_pathway)]
# matrix <- dcast(target, tumor_type~gene_pathway,value.var = "each_tumor_enrich_scores")
# matrix <- setDF(matrix)
# rownames(matrix) <- matrix$tumor_type
# matrix <- matrix[,-1]
# matrix <- t(matrix)
# library(ComplexHeatmap)
# library(circlize)
# col_fun = colorRamp2(c( -5,0,5), c("blue", "white", "red"))
# 
# heatmap_object <- Heatmap(matrix,col = col_fun,
#         cluster_rows = FALSE, cluster_columns = FALSE)
# 
# heatmap_object <- draw(heatmap_object, heatmap_legend_side  = "bottom");
# 
# dev.off()
### reshape heatmap end

#### Breeds within each tumor type 
# need at least 10 dogs,
# need at least two certain dogs have that mutation (gene or variants)

###
base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"
total_sum <- fread(paste(base_dir,"02_25","no_UCL_pathway_summary_02_25.txt",sep = seperator))
## Append Breeds info ##
breed <- match_vector_table(total_sum$sample_names, "final_breed_label", whole_wes_clean_breed_table)
total_sum$Breeds <- breed
total_sum <- total_sum[!Breeds %in% c(NA,"No breed provided and unable to do the breed-prediction"),]

number_breeds_cutoff <- 10
number_sample_mut_cutoff <- 2
### breedwide variants ##
tumor_type <- unique(total_sum$Subtype)
total_tumor_type_summary <- NULL
for (index in 1:length(tumor_type)) {
  each_tumor <- tumor_type[index]
  each_tumor_sum <- NULL
  print(paste( "Processing the ",index," tumor with total tumor",length(tumor_type),sep = " "))
  # each tumor type
  each_tumor_info <- total_sum[Subtype == each_tumor,]
  if (nrow(each_tumor_info) > 0) {
    ## check breeds number
    each_tumor_all_breeds <-
      unique(each_tumor_info[, .(sample_names, Breeds)])
    each_tumor_breed_number <-
      data.table(table(each_tumor_all_breeds$Breeds))
    candidate_breeds <-
      sort(each_tumor_breed_number[N >= number_breeds_cutoff,]$V1)
    candidate_breeds <- candidate_breeds[!candidate_breeds %in% c("Mixed")]
    
    total_breed_total_pathway <- NULL
    
    each_tumor_candidate_breed_info <- each_tumor_info[Breeds %in% candidate_breeds,]
    
    each_breed_each_pathway <- NULL
    for (index in 1:length(path_col)) {
      #print(paste("Processing the ",index," gene mutation with total gene mutation",length(candidate_breeds_uniq_gene_mut),sep = " "))
      each_pathway = path_col[index]
      total_breed_each_pathway <- NULL
      #each_gene_mut = "AKT1"
      #each_tumor_breed_gene <- each_tumor[Breeds ==each_breed]$gene_name
      #each_tumor_breed_uniq_gene <- unique(each_tumor_breed_gene)
      #### check gene mut for each breed ( at least two dogs)
      for (candidate_breed in candidate_breeds){
        total_target_breed_number <- nrow(unique(each_tumor_candidate_breed_info[Breeds==candidate_breed,.(sample_names)]))
        total_other_breeds_number <- nrow(unique(each_tumor_info[Breeds!=candidate_breed,.(sample_names)]))
        target_breed_with <-sum(each_tumor_candidate_breed_info[Breeds==candidate_breed,][[each_pathway]])
        target_breed_without <- total_target_breed_number-target_breed_with
        
        other_breeds_with <- sum(each_tumor_info[Breeds!=candidate_breed,][[each_pathway]])
        other_breeds_without <- total_other_breeds_number-other_breeds_with
        testor <- rbind(c(target_breed_with,target_breed_without),
                        c(other_breeds_with,other_breeds_without ))
        
        fisher_test <- fisher.test(testor)
        p_value <- fisher_test[["p.value"]];
        odds_ratio <- fisher_test[["estimate"]];
        
        enr_score <- min(5, (log10(p_value) * -1));
        
        # if (target_breed_with <number_sample_mut_cutoff){
        #   enr_score <- NA
        # 
        # }
        # else{
          
          if(odds_ratio < 1) {
            enr_score <- 0-enr_score;
          }
          
        # }
        # enrich_p_value <- fisher.test(testor, alternative = "greater")$p.value
        # deplete_p_value <- fisher.test(testor, alternative = "less")$p.value
        # enrich_p_value <- fisher.test(testor, alternative = "greater")$p.value
        # deplete_p_value <- fisher.test(testor, alternative = "less")$p.value
        
        each_breed_each_pathway <- data.table(pathway = each_pathway,
                                              target_breed_with = target_breed_with,
                                              target_breed_without= target_breed_without,
                                              other_breeds_with =other_breeds_with, 
                                              other_breeds_without = other_breeds_without,
                                              enr_score = enr_score,
                                              fisher_pvalue = p_value,
                                              odds_ratio = odds_ratio,
                                              # enrich_p_value = enrich_p_value,
                                              # deplete_p_value=deplete_p_value,
                                              breeds = candidate_breed,
                                              subtype = each_tumor)
        total_breed_each_pathway <- rbindlist(list(total_breed_each_pathway, each_breed_each_pathway))
      }
      
      #if (any(total_breed_each_pathway$target_breed_with >= number_sample_mut_cutoff) ) {
      total_breed_total_pathway <- rbindlist(list(total_breed_total_pathway,total_breed_each_pathway))
      
      #}
      
    }
    
    each_tumor_sum <- rbindlist(list(each_tumor_sum,total_breed_total_pathway))
    
  }
  total_tumor_type_summary <- rbindlist(list(total_tumor_type_summary,each_tumor_sum))
}


#meet_cut_off <- total_tumor_type_summary[target_breed_with >number_sample_mut_cutoff,]

## within each tumor type, do p adjustment for each breed , no need to do BH adjustment because we use enrich scores
# 
# Total_tumor_info <- NULL
# for (each_tumor_type in tumor_type){
#   each_tumor_breed_info <- NULL
#   #each_tumor_type <- "MT"
#   candidate_breed_each_tumor_type <- unique(unlist(total_tumor_type_summary[subtype ==each_tumor_type,.(breeds)]$breeds))
#   for ( each_breed in candidate_breed_each_tumor_type){
#     each_breed_pvalue_for_each_tumor <- total_tumor_type_summary[subtype == each_tumor_type & breeds == each_breed]
#     each_breed_pvalue_for_each_tumor <- each_breed_pvalue_for_each_tumor[order(enrich_p_value)]
#     each_breed_pvalue_for_each_tumor$enrich_BH_pvalue = p.adjust(each_breed_pvalue_for_each_tumor$enrich_p_value, method = "BH")
#     
#     each_breed_pvalue_for_each_tumor <- each_breed_pvalue_for_each_tumor[order(deplete_p_value)]
#     each_breed_pvalue_for_each_tumor$deplete_BH_pvalue = p.adjust(each_breed_pvalue_for_each_tumor$deplete_p_value, method = "BH")
#     
#     each_tumor_breed_info <- rbindlist(list(each_tumor_breed_info, each_breed_pvalue_for_each_tumor))
#   }
#   Total_tumor_info <- rbindlist(list(Total_tumor_info,each_tumor_breed_info))
# }
# target_col <- path_col[-c(2,3,7,8)]

### Reshape the info so that fit the matrix format
target_info <- total_tumor_type_summary[subtype %in% c("MT","TCL","OSA")]
target_info$enr_score
# matrix <- dcast(target_info, breeds+subtype~pathway,value.var = "enr_score")
# matrix <- setDF(matrix)
# rownames(matrix) <- matrix$tumor_type
# matrix <- matrix[,-1]
# matrix <- t(matrix)


fwrite(target_info, file = paste(base_dir,"02_25","with_Pvalue_target_Breeds_sig_pan_tumor_02_25.txt",sep = seperator),
       col.names = T, row.names = F, quote = F,sep = "\t",
       na = "NA")


## Golden retriever compare across different tumor ## 
## golden retriever
# "OSA", "HSA", "BCL", "TCL"

base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"
total_sum <- fread(paste(base_dir,"02_25","no_UCL_pathway_summary_02_25.txt",sep = seperator))
breed <- match_vector_table(total_sum$sample_names, "final_breed_label", whole_wes_clean_breed_table)
total_sum$Breeds <- breed

number_breeds_cutoff <- 10
number_sample_mut_cutoff <- 2
target_breed <- "Golden Retriever"
target_tumor <- c("OSA","BCL","TCL","HSA")

tumor_sum <- NULL
for (tumor in target_tumor){
  target_tumor_target_breed_info <- total_sum[Subtype == tumor & Breeds == target_breed]
  other_tumor_target_breed_info <- total_sum[Subtype != tumor & Breeds == target_breed]
  total_target_tumor_breed_number <- nrow(target_tumor_target_breed_info)
  total_othter_tumor_breed_number <- nrow(total_sum[Subtype != tumor & Breeds == target_breed])
  path_way_sum <- NULL
  
  for (index in 1:length(path_col)){
    each_pathway <- path_col[index]
    target_tumor_target_pathway_with <- sum(target_tumor_target_breed_info[[each_pathway]])
    target_tumor_target_pathway_without <- total_target_tumor_breed_number-target_tumor_target_pathway_with
    other_tumor_target_pathway_with <- sum(other_tumor_target_breed_info[[each_pathway]])
    other_tumor_target_pathway_without <- total_othter_tumor_breed_number - other_tumor_target_pathway_with
    testor <- rbind(c(target_tumor_target_pathway_with,target_tumor_target_pathway_without),
                    c(other_tumor_target_pathway_with,other_tumor_target_pathway_without ))
    
    fisher_test <- fisher.test(testor)
    p_value <- fisher_test[["p.value"]];
    odds_ratio <- fisher_test[["estimate"]];
    
    enr_score <- min(5, (log10(p_value) * -1));
    
    # if (target_tumor_target_pathway_with <number_sample_mut_cutoff){
    #   enr_score <- NA
    #   
    # } 
    # else{
      
      if(odds_ratio < 1) {
        enr_score <- 0-enr_score;
      }
    # }
    each_pathway <- data.table(pathway = each_pathway,
                               tumor_with = target_tumor_target_pathway_with,
                               tumor_without= target_tumor_target_pathway_without,
                               othertumors_with =other_tumor_target_pathway_with, 
                               other_tumors_without = other_tumor_target_pathway_without,
                               enr_score = enr_score,
                               fisher_pvalue = p_value,
                               odds_ratio = odds_ratio,
                               # enrich_p_value = enrich_p_value,
                               # deplete_p_value=deplete_p_value,
                               breeds = target_breed,
                               subtype = tumor)
    path_way_sum <- rbindlist(list(path_way_sum,each_pathway))
    
  }
  tumor_sum <- rbindlist(list(tumor_sum, path_way_sum))
}

fwrite(tumor_sum, file=paste(base_dir,"02_25","GoldenWith_pValue_pathway_02_25.txt",sep=seperator),
       col.names = T, row.names = F, quote = F, sep="\t",eol="\n",na = "NA")

