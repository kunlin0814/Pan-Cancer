library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
seperator <- "/"

#retro_gene_list <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Retro_gene_finding/RetroGeneList/new_retro_gene_list_CanFam3.1.99gtf.txt",
#                         header = F)
whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt")
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt")       
    
# dataset <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Burair_pan_scripts/breed_prediction_test/Pan-Cancer-Breed_prediction/seperate_dis_val/breed_prediction_metadata.txt")

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

amp_delete <- fread(paste(base_dir,"CNV_Taifang_total_amp_delete_no_pseudo_subtype.txt",sep = seperator),
                    header = T)

amp_delete <- amp_delete[!sample_names %in% exclude]
amp_delete <- amp_delete[!grepl("ENSCAFG",amp_delete[,.(gene_name)]$gene_name,ignore.case = T)]
# 
# ## append the column
# breed <- match_vector_table(amp_delete$sample_names, "Breed_info",
#                             table=whole_wes_clean_breed_table, string_value = T)
# amp_delete$breeds <- breed
# subtype <- match_vector_table(amp_delete$sample_names, "DiseaseAcronym2",
#                               table=whole_wes_clean_breed_table, string_value = T)
# amp_delete$subtype <- subtype
# 
# symbol <- match_vector_table(amp_delete$sample_names, "Symbol",
#                              table=whole_wes_clean_breed_table, string_value = T)
# amp_delete$symbol <- symbol
# ## append end
# 
# amp_delete <- amp_delete[,gene_mutation:=paste(gene_name,mut_type,sep = "_")]
unique(amp_delete$sample_names)

# fwrite(amp_delete,file = paste(base_dir,"02_11","CNV_Taifang_total_amp_delete_no_pseudo_subtype.txt",sep = seperator),
#        col.names = T, row.names = F, quote = F, sep = "\t",
#        eol = "\n")

#dataset <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Burair_pan_scripts/breed_prediction_test/Pan-Cancer-Breed_prediction/seperate_dis_val/breed_prediction_metadata.txt")

## append column and write the output
## 
# # colnames(dataset)[2] <- "sample_names"
# # taifang <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp/Taifang/0206/target_col_nonsym_samplewide_p_value_VAF_Mutect_orientBias3_02_01.txt")
# # subtype <- match_vector_table(taifang$sample_names, "DiseaseAcronym2",
# #                               table=whole_wes_clean_breed_table, string_value = T)
# #
# #
# # datatype <- match_vector_table(taifang$sample_names, "Dataset",
# #                               table=dataset, string_value = T)
# # taifang$DataType <- datatype
# 
# fwrite(taifang,file = paste(base_dir,"02_05","target_col_nonsym_samplewide_p_value_VAF_Mutect_orientBias3_02_05.txt",sep = seperator),
#        col.names = T, row.names = F, quote = F, sep = "\t",
#        eol = "\n")

# table_total_sample <- amp_delete$sample_names
# subtype <- sapply(table_total_sample,FUN = match_table, column="DiseaseAcronym2",table=whole_wes_clean_breed_table)
# amp_delete$subtype <- subtype

#breed <- sapply(amp_delete$sample_names,FUN = match_table, column="Breed_info",table=whole_wes_clean_breed_table)

# breed <- match_vector_table(amp_delete$sample_names, "Breed_info", 
#                             table=whole_wes_clean_breed_table, string_value = T)
# amp_delete$Breeds <- breed
# 
# fwrite(amp_delete,file = paste(base_dir,"02_05","With_subtype_no_pesudo_samples_amp_delete_0205.txt",sep = seperator),
#        col.names = T, row.names = F, quote = F, sep = "\t",
#        eol = "\n")

################ main code ####################
tumor_type <- sort(unique(amp_delete$subtype))
## use which to find index and check how many samples ( can't use unique gene to find)
## 2/2 merge amp and delete
## 2/3 seperate merge amp and delete and treat them as different mut


total_gene_summary <- NULL
for (each_tumor in tumor_type) {
  current <- match(each_tumor, tumor_type,nomatch = 0)
  print(paste("Processing the ", current, "tumors with total numbers of tumor", length(tumor_type),sep = " "))
  each_tumor_info_sum <- NULL
  each_tumor_info <- amp_delete[subtype == each_tumor, ]
  each_tumor_mut_type <- each_tumor_info$gene_mutation
  each_tumor_uniq_mut_type <- unique(each_tumor_mut_type)
  each_tumor_total_sample <-   length(unique(each_tumor_info$sample_names))
  
  tumor_type_sum <- character(length(each_tumor_uniq_mut_type))
  gene_mut_sum <- character(length(each_tumor_uniq_mut_type))
  gene_sum <- character(length(each_tumor_uniq_mut_type))
  numbersamples_withgene_sum <-numeric(length(each_tumor_uniq_mut_type))
  numbersamples_withoutgene_sum <- numeric(length(each_tumor_uniq_mut_type))
  total_others_sum <- numeric(length(each_tumor_uniq_mut_type))
  total_others_without_sum <- numeric(length(each_tumor_uniq_mut_type))
  
  total_numbersamples_withgene <- 0
  total_numbersamples_withoutgene <- 0 
  
  # calculate total x and total y
  for (index in 1:length(each_tumor_uniq_mut_type)) {
    #current <- match(each_gene, each_tumor_gene,nomatch = 0)
    print(paste("Processing the ", index, "gene mutation, with total gene mutations", length(each_tumor_uniq_mut_type),sep = " "))
    each_gene <- each_tumor_uniq_mut_type[index]
    total_others <- 0
    total_others_without <- 0
    sample_loc <- which(each_gene == each_tumor_mut_type)
    #length(unique(each_tumor_info[gene_mutation == each_gene, ]$sample_names))
    numbersamples_withgene <- length(unique(each_tumor_info$sample_names[sample_loc]))
    #length(which(each_gene == each_tumor_gene))
    #length(unique(each_tumor_info[gene_mutation == each_gene, ]$sample_names))
    
    numbersamples_withoutgene <- each_tumor_total_sample- numbersamples_withgene
    outside_target_gene <- each_tumor_uniq_mut_type[!each_tumor_uniq_mut_type %in% each_gene]
    
    total_numbersamples_withgene <- total_numbersamples_withgene+numbersamples_withgene
    total_numbersamples_withoutgene <- total_numbersamples_withoutgene+numbersamples_withoutgene
    
    #gene_sum[index] <- each_gene
    tumor_type_sum[index] <- each_tumor
    gene_mut_sum[index] <- each_gene
    numbersamples_withgene_sum[index] <- numbersamples_withgene
    numbersamples_withoutgene_sum[index] <- numbersamples_withoutgene
    
    
  }
  total_others_sum <- total_numbersamples_withgene -numbersamples_withgene_sum
  total_others_without_sum <- total_numbersamples_withoutgene- numbersamples_withoutgene_sum
  
  
  
  
  each_tumor_info_sum <- data.table( tumor_type = tumor_type_sum,
                                     gene_mut = gene_mut_sum,
                                     numbersamples_withgene = numbersamples_withgene_sum,
                                     numbersamples_withoutgene = numbersamples_withoutgene_sum,
                                     total_others = total_others_sum,
                                     total_others_withoutgene = total_others_without_sum
  )
  
  p_value_sum <- numeric(length(each_tumor_uniq_mut_type))
  for ( i in 1:nrow(each_tumor_info_sum)) {
    
    target <- as.matrix(each_tumor_info_sum[i,.(numbersamples_withgene,numbersamples_withoutgene)])
    others <- as.matrix(each_tumor_info_sum[i,.(total_others,total_others_withoutgene)])
    testor <- rbind(c(target),c(others))
    each_gene_p_value <- fisher.test(testor, alternative = "greater")$p.value
    p_value_sum[i] <-  each_gene_p_value
    
  }
  each_tumor_info_sum$p_value <- p_value_sum 
  total_gene_summary <- rbindlist(list(total_gene_summary, each_tumor_info_sum))
}



#### Breeds within each tumor type 
# need at least 10 dogs,
# need at least two certain dogs have that mutation (gene or variants)
base_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
seperator <- "/"

#retro_gene_list <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Retro_gene_finding/RetroGeneList/new_retro_gene_list_CanFam3.1.99gtf.txt",
#                         header = F)
whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt")
#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt")       

# dataset <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Burair_pan_scripts/breed_prediction_test/Pan-Cancer-Breed_prediction/seperate_dis_val/breed_prediction_metadata.txt")

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

amp_delete <- fread(paste(base_dir,"CNV_Taifang_total_amp_delete_no_pseudo_subtype.txt",sep = seperator),
                    header = T)

amp_delete <- amp_delete[!sample_names %in% exclude]
amp_delete <- amp_delete[!grepl("ENSCAFG",amp_delete[,.(gene_name)]$gene_name,ignore.case = T)]

amp_delete$

number_breeds_cutoff <- 10
number_sample_mut_cutoff <- 2
### breedwide variants ##
tumor_type <- unique(mutect_after_vaf$Subtype)
total_tumor_type_summary <- NULL
for (index in 1:length(tumor_type)) {
  each_tumor <- tumor_type[index]
  each_tumor_sum <- NULL
  print(paste(
    "Processing the ",
    index,
    " tumor with total tumor",
    length(tumor_type),
    sep = " "
  ))
  # each tumor type
  each_tumor_info <-
    mutect_after_vaf[Subtype == each_tumor &
                       !Breeds %in% c("Mixed", NA)]
  if (nrow(each_tumor_info) > 0) {
    ## check breeds number
    each_tumor_all_breeds <-
      unique(each_tumor_info[, .(sample_names, Breeds)])
    each_tumor_breed_number <-
      data.table(table(each_tumor_all_breeds$Breeds))
    candidate_breeds <-
      sort(each_tumor_breed_number[N >= number_breeds_cutoff,]$V1)
    candidate_breeds_uniq_gene_mut <-
      unique(mutect_after_vaf[Breeds %in% candidate_breeds, .(gene_name)]$gene_name)
    all_gene_summary <- NULL
    for (index in 1:length(candidate_breeds_uniq_gene_mut)) {
      #print(paste("Processing the ",index," gene mutation with total gene mutation",length(candidate_breeds_uniq_gene_mut),sep = " "))
      each_gene_mut = candidate_breeds_uniq_gene_mut[index]
      #each_gene_mut = "AKT1"
      #each_tumor_breed_gene <- each_tumor[Breeds ==each_breed]$gene_name
      #each_tumor_breed_uniq_gene <- unique(each_tumor_breed_gene)
      #### check gene mut for each breed ( at least two dogs)
      candidate_breeds_info <-
        unique(each_tumor_info[gene_name == each_gene_mut, .(sample_names, Breeds)])
      number_gene_mut_in_breeds <-
        candidate_breeds_info[, .N, keyby = .(Breeds)]
      if (any(number_gene_mut_in_breeds$N >= number_sample_mut_cutoff) ) {
        missing_breed <-
          setdiff(candidate_breeds, number_gene_mut_in_breeds$Breeds)
        if (length(missing_breed) > 0) {
          for (each_missing in missing_breed) {
            missing_info <- list(Breeds = each_missing, N = 0)
            number_gene_mut_in_breeds <-
              rbindlist(list(number_gene_mut_in_breeds, missing_info))
          }
        }
        number_gene_mut_in_breeds <-
          number_gene_mut_in_breeds[Breeds %in% candidate_breeds]
        total_candidate_dogs <-
          length(unique(each_tumor_info[Breeds %in% candidate_breeds]$sample_names))
        total_candidate_breeds_with <- sum(number_gene_mut_in_breeds$N)
        total_candidate_breeds_without <-
          total_candidate_dogs - total_candidate_breeds_with
        
        all_candidate_breed_each_gene_sum <- NULL
        #each_tumor_sum <- rbindlist(list(each_tumor_sum,number_gene_mut_in_breeds))
        for (each_candidate_breed in candidate_breeds) {
          total_target_sample <-
            nrow(unique(each_tumor_info[Breeds == each_candidate_breed, .(sample_names, Breeds)]))
          target_with <-
            number_gene_mut_in_breeds[Breeds == each_candidate_breed, ]$N
          others_with <- total_candidate_breeds_with - target_with
          target_without <- total_target_sample - target_with
          others_without <- total_candidate_breeds_without - target_without
          testor <-
            rbind(c(target_with, target_without),
                  c(others_with, others_without))
          
          each_breed_p_value <-
            fisher.test(testor, alternative = "greater")$p.value
          
          each_breed_sum <-
            data.table(
              breeds = each_candidate_breed,
              gene_mut = each_gene_mut,
              tumor_type = each_tumor,
              target_breeds_with = target_with,
              target_breeds_without = target_without,
              others_breeds_with = others_with,
              others_breeds_without = others_without,
              p_value = each_breed_p_value
            )
          
          all_candidate_breed_each_gene_sum <-
            rbindlist(list(all_candidate_breed_each_gene_sum, each_breed_sum))
        }
        all_candidate_breed_each_gene_sum <- setDT(all_candidate_breed_each_gene_sum)
        #all_candidate_breed_each_gene_sum <- all_candidate_breed_each_gene_sum[order(p_value)]
        #all_candidate_breed_each_gene_sum$BH_pvalue = p.adjust(all_candidate_breed_each_gene_sum$p_value, method = "BH")
        all_gene_summary <-
          rbindlist(list(all_gene_summary, all_candidate_breed_each_gene_sum))
      }
    }
    each_tumor_sum <-
      rbindlist(list(each_tumor_sum, all_gene_summary))
  }
  total_tumor_type_summary <-
    rbindlist(list(total_tumor_type_summary, each_tumor_sum))
}

# fwrite(total_tumor_type_summary, file = "C:/Users/abc73/Desktop/Breed_associated_sig_pvalue_02_13.txt",
#        col.names = T, row.names = F, quote = F, sep = "\t")


# total_tumor_type_summary <- fread("C:/Users/abc73/Desktop/Breed_associated_sig_pvalue_02_13.txt")
unique(total_tumor_type_summary$tumor_type)

meet_cut_off <- total_tumor_type_summary[target_breeds_with >number_sample_mut_cutoff,]

## within each tumor type, do p adjustment for each breed

Total_tumor_info <- NULL
for (each_tumor_type in tumor_type){
  each_tumor_breed_info <- NULL
  #each_tumor_type <- "MT"
  candidate_breed_each_tumor_type <- unique(unlist(meet_cut_off[tumor_type ==each_tumor_type,.(breeds)]$breeds))
  for ( each_breed in candidate_breed_each_tumor_type){
    each_breed_pvalue_for_each_tumor <- meet_cut_off[tumor_type == each_tumor_type & breeds == each_breed]
    each_breed_pvalue_for_each_tumor <- each_breed_pvalue_for_each_tumor[order(p_value)]
    each_breed_pvalue_for_each_tumor$BH_pvalue = p.adjust(each_breed_pvalue_for_each_tumor$p_value, method = "BH")
    each_tumor_breed_info <- rbindlist(list(each_tumor_breed_info, each_breed_pvalue_for_each_tumor))
  }
  Total_tumor_info <- rbindlist(list(Total_tumor_info,each_tumor_breed_info))
}

split_col <- as.data.frame(str_split_fixed(Total_tumor_info$gene_mut,"_",2))
colnames(split_col) <- c("gene_name","mut_type")
Total_tumor_info <- cbind(Total_tumor_info,split_col)

fwrite(Total_tumor_info, file = paste(base_dir,"WithBH_amp_delete_breed_significant_Tumor_wide_02_15.txt",sep=seperator),
       col.names = T, row.names = F, quote = F, eol = "\n",na = "NA",
       sep = "\t")


### breed tumor wide end ###




### Analyzed the results, sep gene and mut type

amp_delete <- total_gene_summary

different_type <- as.data.table(str_split_fixed(amp_delete$gene_mut,"_",2))
colnames(different_type) <- c("gene_name","mut_type")
amp_delete <- cbind(amp_delete,different_type)


#### adjust pvalue with BH methods ####


total_summary <- NULL
tumor_type <- sort(unique(amp_delete$tumor_type))

for (each_tumor in tumor_type) {
    each_tumor_sum <- NULL
    each_tumor_info <- amp_delete[tumor_type==each_tumor]
    #gene_name <- sort(unique(each_tumor_info[,.(gene_name)]$gene_name))
    #each_mut_type_tumor <- each_tumor_info[gene_name==each_mut_type,]
    each_mut_type_tumor <- each_tumor_info[order(p_value)]
    each_mut_type_tumor$BH_pvalue <- p.adjust(each_mut_type_tumor$p_value, method = "BH")
    each_tumor_sum <- rbindlist(list(each_tumor_sum,each_mut_type_tumor))
    

  total_summary <- rbindlist(list(total_summary,each_tumor_sum))
}

# fwrite(total_summary, file = paste(base_dir,"02_11", "With_BH_sep_amp_delete_merged_adjust_Pan_cancer_amp_delete_02_11.txt",
#                                         sep = seperator),
#        col.names = T, row.names = F, quote = F,sep = "\t",
#        eol = "\n")

### Analyzed the final results ###
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)
base_dir <- "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Gene_amp"
seperator <- "/"

original_amp_delete <- fread(paste(base_dir, "With_subtype_no_pesudo_samples_amp_delete_0202.txt", sep = seperator))

p_value_amp_delete <- fread(paste(base_dir,"02_03", "With_BH_sep_amp_delete_merged_adjust_Pan_cancer_amp_delete_02_03.txt", sep = seperator))                                    
merge_gene_amp_delete <- fread(paste(base_dir,"02_02", "With_BH_final_Pan_cancer_amp_delete_02_02.txt", sep = seperator))                                    

# check the result
top10_sep <- p_value_amp_delete[,head(.SD,10),keyby = tumor_type] 
top10_merge <- merge_gene_amp_delete[,head(.SD,10),keyby = tumor_type]
MT <- p_value_amp_delete[tumor_type=="MT" ]
a <- MT[grepl("CDKN",p_value_amp_delete[tumor_type=="MT",.(gene_name)][["gene_name"]], ignore.case = T)]

sig_amp_del <- p_value_amp_delete[BH_pvalue < 0.05]

# mut_type <- amp_delete[ ,.N, keyby = .(tumor_type,mut_type)]



