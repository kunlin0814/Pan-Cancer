library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")

base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt") 

mutect_after_vaf <- fread(paste(base_dir,"NonSyn_Burair_filtering3_WithBreeds_Subtypes_QCpass_mutect_after_vaf_02_11.txt",
                                sep =seperator))
amp_delete <- fread(paste(base_dir,"CNV_Taifang_total_amp_delete_no_pseudo_subtype.txt",sep = seperator),header = T)
SNV <- unique(mutect_after_vaf[,c("sample_names","gene_name","status","Subtype"), with =F])
CNV <- unique(amp_delete[,c("sample_names","gene_name","mut_type","subtype"),with =F])
colnames(CNV)<- c("sample_names","gene_name","status","Subtype")
total_mut <- rbindlist(list(SNV,CNV))
total_sample <- unique(total_mut$sample_names)
pathway <- fread(paste(base_dir,"all_pathway.txt",sep = seperator), na.strings = "")
path_col <- colnames(pathway)

# need to put sample wide and tumor wide

Other_PI3K_AKT <- numeric(length(total_sample))
MAPK <- numeric(length(total_sample))
RTK <- numeric(length(total_sample))
Other_Cell_cycle <- numeric(length(total_sample))
p53 <- numeric(length(total_sample))
BCR <- numeric(length(total_sample))
Chromatin_remodeler <- numeric(length(total_sample))
Other_p53 <- numeric(length(total_sample))
CDKN2A <- numeric(length(total_sample))
TCR <- numeric(length(total_sample))
MDM2 <- numeric(length(total_sample))
PIK3CA <- numeric(length(total_sample))
sample_names <- character(length(total_sample))
total_sum <- list(sample_names=sample_names, 
                  Other_PI3K_AKT= Other_PI3K_AKT,
                  MAPK = MAPK,
                  RTK = RTK,
                  Other_Cell_cycle = Other_Cell_cycle,
                  p53 = p53,
                  BCR = BCR,
                  Chromatin_remodeler = Chromatin_remodeler,
                  Other_p53 = Other_p53,
                  CDKN2A = CDKN2A,
                  TCR = TCR,
                  MDM2 = MDM2,
                  PIK3CA = PIK3CA)

#each_sample_gene <- total_mut[sample_names=="004"]$gene_name
for (index in 1:length(total_sample)){
  each_sample = total_sample[index]
  total_sum[["sample_names"]][index] = each_sample
  each_sample_gene <- total_mut[sample_names==each_sample]$gene_name
  
  for ( col_index in 1:length(path_col)){
    col_name <- path_col[col_index]
    each_path_way_gene <- pathway[, col_name, with = F][[col_name]]
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

total_tumor_p_value_sum <- NULL
no_UCL <- total_sum[Subtype!="UCL",]

fwrite(no_UCL,file = paste(base_dir,"no_UCL_pathway_summary.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t")

total_tumor_type <- unique(no_UCL$Subtype)
for (index in 1:length(total_tumor_type)){
  
  each_tumor <- total_tumor_type[index]
  total_target_sample_number <- length(unique(no_UCL[Subtype==each_tumor]$sample_name))
  total_other_sample_number <- length(unique(no_UCL[Subtype!=each_tumor]$sample_name))
  each_tumor_target <- numeric(length(path_col))
  each_tumor_target_without <- numeric(length(path_col))
  each_tumor_other_with <- numeric(length(path_col))
  each_tumor_other_without <- numeric(length(path_col))
  each_tumor_enrich_p_value <- numeric(length(path_col))
  each_tumor_deplete_p_value <- numeric(length(path_col))
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
    
    enrich_p_value <- fisher.test(testor, alternative = "greater")$p.value
    deplete_p_value <- fisher.test(testor, alternative = "less")$p.value
    
    each_tumor_enrich_p_value[col_index] <- enrich_p_value
    each_tumor_deplete_p_value[col_index] <- deplete_p_value
    #print(each_col)
  }
  each_tumor_sum <- data.table(tumor_type = each_tumor_type,
                               gene_pathway = each_tumor_col,
                               each_gene_with = each_tumor_target,
                               each_gene_without = each_tumor_target_without,
                               each_gene_other_with = each_tumor_other_with,
                               each_gene_other_without = each_tumor_other_without,
                               each_gene_enrich_pvalue = each_tumor_enrich_p_value,
                               each_gene_deplete_pvalue= each_tumor_deplete_p_value
  )
  
  each_tumor_sum <- each_tumor_sum[order(each_gene_enrich_pvalue)]
  each_tumor_sum$Enrich_BH_pvalue = p.adjust(each_tumor_sum$each_gene_enrich_pvalue, method = "BH")
  each_tumor_sum <- each_tumor_sum[order(each_gene_deplete_pvalue)]
  each_tumor_sum$Deplete_BH_pvalue = p.adjust(each_tumor_sum$each_gene_deplete_pvalue, method = "BH")

  total_tumor_p_value_sum <- rbind(total_tumor_p_value_sum, each_tumor_sum)
  
}

# fwrite(total_tumor_p_value_sum, file=paste(base_dir,"02_11","Pathway","pan-tumor_pathway_02_11.txt",sep=seperator),
#        col.names = T, row.names = F, quote = F, sep="\t",eol="\n")



#### reshape and heatmap ###
# 
# target <- total_tumor_p_value_sum[,.(tumor_type,Enrich_BH_pvalue,gene_pathway)]
# matrix <- dcast(target, tumor_type~gene_pathway,value.var = "Enrich_BH_pvalue")
# matrix <- setDF(matrix)
# rownames(matrix) <- matrix$tumor_type
# matrix <- matrix[,-1]
# matrix <- t(matrix)
# library(ComplexHeatmap)
# 
# 
# 
# col_fun = colorRamp2(c( 0,0.5,1), c("green", "white", "red"))
# 
# Heatmap(matrix,col = col_fun,
#         cluster_rows = FALSE, cluster_columns = FALSE)




### reshape heatmap end

#### Breeds within each tumor type 
# need at least 10 dogs,
# need at least two certain dogs have that mutation (gene or variants)

## Append Breeds info ##

# breed <- match_vector_table(total_mut$sample_names, "Breed_info", whole_wes_clean_breed_table)
# total_mut$Breeds <-breed 

breed <- match_vector_table(total_sum$sample_names, "Breed_info", whole_wes_clean_breed_table)
total_sum$Breeds <- breed
#total_mut <- total_mut[Subtype!="UCL",]

###
base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"
total_sum <- fread(paste(base_dir,"no_UCL_pathway_summary.txt",sep = seperator))
breed <- match_vector_table(total_sum$sample_names, "Breed_info", whole_wes_clean_breed_table)
total_sum$Breeds <- breed

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
  each_tumor_info <- total_sum[Subtype == each_tumor & !Breeds %in% c("Mixed", NA)]
  if (nrow(each_tumor_info) > 0) {
    ## check breeds number
    each_tumor_all_breeds <-
      unique(each_tumor_info[, .(sample_names, Breeds)])
    each_tumor_breed_number <-
      data.table(table(each_tumor_all_breeds$Breeds))
    candidate_breeds <-
      sort(each_tumor_breed_number[N >= number_breeds_cutoff,]$V1)
    total_breed_total_pathway <- NULL
    
    each_tumor_info <- each_tumor_info[Breeds %in% candidate_breeds,]
    
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
        total_target_breed_number <- nrow(unique(each_tumor_info[Breeds==candidate_breed,.(sample_names)]))
        total_other_breeds_number <- nrow(unique(each_tumor_info[Breeds!=candidate_breed,.(sample_names)]))
        
        target_breed_with <-sum(each_tumor_info[Breeds==candidate_breed,][[each_pathway]])
        target_breed_without <- total_target_breed_number-target_breed_with
        
        other_breeds_with <- sum(each_tumor_info[Breeds!=candidate_breed,][[each_pathway]])
        other_breeds_without <- total_other_breeds_number-other_breeds_with
        testor <- rbind(c(target_breed_with,target_breed_without),
                        c(other_breeds_with,other_breeds_without ))
        
        enrich_p_value <- fisher.test(testor, alternative = "greater")$p.value
        deplete_p_value <- fisher.test(testor, alternative = "less")$p.value
        
        each_breed_each_pathway <- data.table(pathway = each_pathway,
                                              target_breed_with = target_breed_with,
                                              target_breed_without= target_breed_without,
                                              other_breeds_with =other_breeds_with, 
                                              other_breeds_without = other_breeds_without,
                                              enrich_p_value = enrich_p_value,
                                              deplete_p_value=deplete_p_value,
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

tumor_type <- unique(total_tumor_type_summary$subtype)

## within each tumor type, do p adjustment for each breed

Total_tumor_info <- NULL
for (each_tumor_type in tumor_type){
  each_tumor_breed_info <- NULL
  #each_tumor_type <- "MT"
  candidate_breed_each_tumor_type <- unique(unlist(total_tumor_type_summary[subtype ==each_tumor_type,.(breeds)]$breeds))
  for ( each_breed in candidate_breed_each_tumor_type){
    each_breed_pvalue_for_each_tumor <- total_tumor_type_summary[subtype == each_tumor_type & breeds == each_breed]
    each_breed_pvalue_for_each_tumor <- each_breed_pvalue_for_each_tumor[order(enrich_p_value)]
    each_breed_pvalue_for_each_tumor$enrich_BH_pvalue = p.adjust(each_breed_pvalue_for_each_tumor$enrich_p_value, method = "BH")
    
    each_breed_pvalue_for_each_tumor <- each_breed_pvalue_for_each_tumor[order(deplete_p_value)]
    each_breed_pvalue_for_each_tumor$deplete_BH_pvalue = p.adjust(each_breed_pvalue_for_each_tumor$deplete_p_value, method = "BH")
    
    each_tumor_breed_info <- rbindlist(list(each_tumor_breed_info, each_breed_pvalue_for_each_tumor))
  }
  Total_tumor_info <- rbindlist(list(Total_tumor_info,each_tumor_breed_info))
}
path_col
target_col <- path_col[-c(2,3,7,8)]

target_info <- Total_tumor_info[subtype %in% c("MT","BCL","TCL","OSA") & pathway %in% target_col,]
fwrite(Total_tumor_info, file = paste(base_dir,"02_11","Pathway","total_Breeds_sig_pan_tumor.txt",sep = seperator),
       col.names = T, row.names = F, quote = F,sep = "\t")


