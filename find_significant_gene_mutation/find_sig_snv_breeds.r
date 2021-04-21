#### Breeds within each tumor type 
# need at least 10 dogs,
# need at least two certain dogs have that mutation (gene or variants)
base_dir <- 
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1/"
seperator <- "/"
whole_wes_clean_breed_table <- fread(#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt") 
  "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt")

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))


SNV <- fread(paste(base_dir,"Final_Total_withGene_Burair_Filtering3_VAF_Mutect_orientBiasModified_04_02.txt.gz",sep =seperator)
             ,fill = T,na.strings="")

SNV <- SNV[status!= "synonymous",]
SNV <- SNV[!sample_names %in% exclude]
SNV <- SNV[tumor_type!="UCL"]
SNV$Subtype <- match_vector_table(SNV$sample_names,column="DiseaseAcronym2",table=whole_wes_clean_breed_table)

SNV$Breeds  <- match_vector_table(SNV$sample_names,column="final_breed_label",table=whole_wes_clean_breed_table)

#total_breeds <- unique( SNV$Breeds)
# clean_breeds <- na.omit(total_breeds)
indel_file <- fread(paste(base_dir,"total_CDS_indel_info_withGene_04_08.txt",sep =seperator))
colnames(indel_file) <- c('chrom','pos','ref','alt','gene_name','ensembl_id','status','sample_names')
indel_file <- indel_file[gene_name!="-" & status!="nonframeshift " & ! sample_names %in% exclude,]
indel_file$Subtype <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table)
indel_file$Breeds <- match_vector_table(indel_file$sample_names,"final_breed_label",whole_wes_clean_breed_table)


#setcolorder(indel_file,c("sample_names","gene_name","emsembl_id","status"))
#indel_file <- indel_file[,emsembl_id:=NULL]


final_indel <- indel_file[, .(sample_names,gene_name,status,Subtype,Breeds)]
final_SNV <- SNV[,.(sample_names,gene_name,status,Subtype,Breeds)]

snv_indel <- rbindlist(list(final_SNV,final_indel))

snv_indel <- snv_indel[Subtype!="UCL" & !sample_names %in% exclude & Breeds!="No breed provided and unable to do the breed-prediction",]
snv_indel <- na.omit(snv_indel)
#### Breeds within each tumor type 
# need at least 10 dogs,
# need at least two certain dogs have that mutation (gene or variants)

number_breeds_cutoff <- 10
number_sample_mut_cutoff <- 2
### breedwide variants ##
tumor_type <- unique(snv_indel$Subtype)
total_tumor_type_summary <- NULL
for (index in 1:length(tumor_type)) {
  each_tumor <- tumor_type[index]
  #each_tumor <- "OM"
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
    snv_indel[Subtype == each_tumor ]
  
  if (nrow(each_tumor_info) > 0) {
    ## check breeds number
    each_tumor_all_breeds <-
      unique(each_tumor_info[, .(sample_names, Breeds)])
    each_tumor_breed_number <-
      data.table(table(each_tumor_all_breeds$Breeds))
    candidate_breeds <-
      sort(each_tumor_breed_number[N >= number_breeds_cutoff,]$V1)
    candidate_breeds <- candidate_breeds[!candidate_breeds %in% c("Mixed")]
    all_gene_summary <- NULL
    if (length(candidate_breeds)>=2){
      candidate_breeds_uniq_gene_mut <-
        unique(snv_indel[Breeds %in% candidate_breeds, .(gene_name)]$gene_name)
      
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
          ####### here total candidate dogs become total dogs #####
          total_candidate_dogs <-
            length(unique(each_tumor_info$sample_names))
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

meet_cut_off <- total_tumor_type_summary[target_breeds_with >number_sample_mut_cutoff,]

## within each tumor type, do p adjustment for each breed
library(data.table)
library(tidyverse)
library(readxl)
source(#"C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
  "/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- 
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1/"
seperator <- "/"

whole_wes_clean_breed_table <- fread(#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt") 
  "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt")

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))
output_dir <- paste(base_dir,"04_07",sep = seperator)

#### Breeds within each tumor type 
# need at least 10 dogs,
# need at least two certain dogs have that mutation (gene or variants)
base_dir <- 
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1/"
seperator <- "/"
whole_wes_clean_breed_table <- fread(#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt") 
  "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt")

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))


SNV <- fread(paste(base_dir,"Final_Total_withGene_Burair_Filtering3_VAF_Mutect_orientBiasModified_04_02.txt.gz",sep =seperator)
             ,fill = T,na.strings="")

SNV <- SNV[status!= "synonymous",]
SNV <- SNV[!sample_names %in% exclude]
SNV <- SNV[tumor_type!="UCL"]
SNV$Subtype <- match_vector_table(SNV$sample_names,column="DiseaseAcronym2",table=whole_wes_clean_breed_table)

SNV$Breeds  <- match_vector_table(SNV$sample_names,column="final_breed_label",table=whole_wes_clean_breed_table)

#total_breeds <- unique( SNV$Breeds)
# clean_breeds <- na.omit(total_breeds)
indel_file <- fread(paste(base_dir,"total_CDS_indel_info_withGene_04_08.txt",sep =seperator))
colnames(indel_file) <- c('chrom','pos','ref','alt','gene_name','ensembl_id','status','sample_names')
indel_file <- indel_file[gene_name!="-" & status!="nonframeshift " & ! sample_names %in% exclude,]
indel_file$Subtype <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table)
indel_file$Breeds <- match_vector_table(indel_file$sample_names,"final_breed_label",whole_wes_clean_breed_table)
final_indel <- indel_file[, .(sample_names,gene_name,status,Subtype,Breeds)]
final_SNV <- SNV[,.(sample_names,gene_name,status,Subtype,Breeds)]

snv_indel <- rbindlist(list(final_SNV,final_indel))

snv_indel <- snv_indel[Subtype!="UCL" & !sample_names %in% exclude & Breeds!="No breed provided and unable to do the breed-prediction",]
snv_indel <- na.omit(snv_indel)


#### Breeds within each tumor type 
# need at least 10 dogs,
# need at least two certain dogs have that mutation (gene or variants)

number_breeds_cutoff <- 10
number_sample_mut_cutoff <- 2
### breedwide variants ##
tumor_type <- unique(snv_indel$Subtype)
total_tumor_type_summary <- NULL
for (index in 1:length(tumor_type)) {
  each_tumor <- tumor_type[index]
  #each_tumor <- "OM"
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
    snv_indel[Subtype == each_tumor ]
  
  if (nrow(each_tumor_info) > 0) {
    ## check breeds number
    each_tumor_all_breeds <-
      unique(each_tumor_info[, .(sample_names, Breeds)])
    each_tumor_breed_number <-
      data.table(table(each_tumor_all_breeds$Breeds))
    candidate_breeds <-
      sort(each_tumor_breed_number[N >= number_breeds_cutoff,]$V1)
    candidate_breeds <- candidate_breeds[!candidate_breeds %in% c("Mixed")]
    all_gene_summary <- NULL
    if (length(candidate_breeds)>=2){
      candidate_breeds_uniq_gene_mut <-
        unique(snv_indel[Breeds %in% candidate_breeds, .(gene_name)]$gene_name)
      
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
          ####### here total candidate dogs become total dogs #####
          total_candidate_dogs <-
            length(unique(each_tumor_info$sample_names))
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

meet_cut_off <- total_tumor_type_summary[target_breeds_with >number_sample_mut_cutoff,]

## within each tumor type, do p adjustment for each breed

Total_tumor_info <- NULL
for (each_tumor_type in tumor_type){
  each_tumor_final_breed_label_label <- NULL
  #each_tumor_type <- "MT"
  candidate_breed_each_tumor_type <- unique(unlist(meet_cut_off[tumor_type ==each_tumor_type,.(breeds)]$breeds))
  for ( each_breed in candidate_breed_each_tumor_type){
    each_breed_pvalue_for_each_tumor <- meet_cut_off[tumor_type == each_tumor_type & breeds == each_breed]
    each_breed_pvalue_for_each_tumor <- each_breed_pvalue_for_each_tumor[order(p_value)]
    each_breed_pvalue_for_each_tumor$BH_pvalue = p.adjust(each_breed_pvalue_for_each_tumor$p_value, method = "BH")
    each_tumor_final_breed_label_label <- rbindlist(list(each_tumor_final_breed_label_label, each_breed_pvalue_for_each_tumor))
  }
  Total_tumor_info <- rbindlist(list(Total_tumor_info,each_tumor_final_breed_label_label))
}

fwrite(Total_tumor_info, file = paste(output_dir,"final_breed_label_sig_WithBH_Tumor_wide_04_07.txt",sep=seperator),
       col.names = T, row.names = F, quote = F, eol = "\n", na = "NA",
       sep = "\t")