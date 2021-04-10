library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
  #"/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- 
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1/"
seperator <- "/"

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt") 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_03_30.txt")

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))
output_dir <- paste(base_dir,"04_07",sep = seperator)

whole_wes_clean_breed_table$final_breed_label
#### Tumor wide analysis ####
### assume all indel are significant
indel_file <- fread(paste(base_dir,"total_CDS_indel_info_withGene_04_08.txt",sep =seperator))
colnames(indel_file) <- c('chrom','pos','ref','alt','gene_name','ensembl_id','status','sample_names')
indel_file <- indel_file[!sample_names %in% exclude,]
#indel_file <- indel_file[gene_name!="-" & status!="nonframeshift " & ! sample_names %in% exclude,]
indel_file$Subtype <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table)
indel_file$Breeds <- match_vector_table(indel_file$sample_names,"final_breed_label",whole_wes_clean_breed_table)
# fwrite(indel_file,
#        file = paste(output_dir,"ExcludeFailQC_CDS_indel_info_withGene_04_08.txt",sep = seperator)
#        ,col.names = T,row.names = F,
#        quote = F,
#        eol = "\n",
#        sep ="\t")

#indel_file$chrom_loc <- paste(indel_file$chrom,indel_file$pos,sep = "_")
# 
# sig_indel_variants <- indel_file[,.(chrom,pos,ref,alt,sample_names,gene_name,ensembl_id,status,
#                                     Breeds)]
# 
# sig_indel_variants$chrom_loc <- paste(sig_indel_variants$chrom,sig_indel_variants$pos,sep = "_")
# sig_indel_variants$BH_pvalue <- 0.01
# sig_indel_variants$Subtype <- match_vector_table(sig_indel_variants$sample_names,'DiseaseAcronym2',whole_wes_clean_breed_table)
# 
# sig_variants_sample_wide <- fread(file = paste(output_dir,"final_breed_variant_nonsyn_samplewide_p_value_orient_modify_04_07.txt",sep = seperator))
# sig_variants_sample_wide <- sig_variants_sample_wide[,.(chrom,pos,ref,alt,sample_names,gene_name,ensembl_id,status,
#                                                         Breeds,chrom_loc,BH_pvalue)]
# 
# sig_variants_sample_wide$Subtype <- match_vector_table(sig_variants_sample_wide$sample_names,'DiseaseAcronym2',whole_wes_clean_breed_table)
# 
# sig_variants_sample_wide <- rbind(sig_variants_sample_wide,sig_indel_variants)
# 
# ### combine snv indel for mutect file
# mutect_after_vaf <- fread(paste(base_dir,"Final_Total_withGene_Burair_Filtering3_VAF_Mutect_orientBiasModified_04_02.txt.gz",sep =seperator),
#                           na.strings = "")
# mutect_after_vaf <- mutect_after_vaf[status!= "synonymous",]
# mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude]
# mutect_after_vaf <- mutect_after_vaf[tumor_type!="UCL"]
# mutect_after_vaf$Subtype <- match_vector_table(mutect_after_vaf$sample_names,column = "DiseaseAcronym2", table =whole_wes_clean_breed_table,string_value = T )
# finalbreed <- match_vector_table(mutect_after_vaf$sample_names,column="final_breed_label",table=whole_wes_clean_breed_table)
# mutect_after_vaf$Breeds <- finalbreed
# mutect_after_vaf <- mutect_after_vaf[,chrom_loc:= paste(chrom,pos,sep = "_"),]
# #mutect_after_vaf <- mutect_after_vaf[,.(sample_names,chrom_loc,Subtype,gene_name,ensembl_id)]
# tumor_type <- unique(mutect_after_vaf$Subtype)
# 
# final_mutect_after_vaf <- mutect_after_vaf[,.(sample_names,gene_name,ensembl_id,status,ref,alt,Subtype,chrom_loc,Breeds)]
# final_indel_file <- indel_file[,.(sample_names,gene_name,ensembl_id,status,ref,alt,Subtype,chrom_loc,
#                                   Breeds)]
# final_variant_tumor_wide_SNV_indel <- rbind(final_mutect_after_vaf,final_indel_file)
# ### combine snv indel for genes end ###
# ## need output
# 
# ## Tumorwide variant analysis
# tumor_type <- unique(final_variant_tumor_wide_SNV_indel$Subtype)
# 
# ### Tumorwide variants
# total_variant_summary <- NULL
# for (each_tumor in tumor_type) {
#   current <- match(each_tumor, tumor_type,nomatch = 0)
#   print(paste("Processing the ", current, "tumors with total numbers of tumor", length(tumor_type),sep = " "))
#   each_tumor_info_sum <- NULL
#   each_tumor_info <- final_variant_tumor_wide_SNV_indel[Subtype == each_tumor, ]
#   each_tumor_variants <- each_tumor_info$chrom_loc
#   
#   # the sig variants come from samplewide vaf comparision 
#   each_tumor_uniq_sig_variants <- unique(sig_variants_sample_wide[BH_pvalue < 0.05 & Subtype==each_tumor, .(chrom_loc)]$chrom_loc)
#   
#   ## vectorized the value I want to append later
#   tumor_type_sum <- character(length(each_tumor_uniq_sig_variants))
#   gene_mut_sum <- character(length(each_tumor_uniq_sig_variants))
#   numbersamples_withgene_sum <-numeric(length(each_tumor_uniq_sig_variants))
#   numbersamples_withoutgene_sum <- numeric(length(each_tumor_uniq_sig_variants))
#   total_others_sum <- numeric(length(each_tumor_uniq_sig_variants))
#   total_others_without_sum <- numeric(length(each_tumor_uniq_sig_variants))
#   each_tumor_total_sample <-   length(unique(each_tumor_info$sample_names))
#   ensmbl_id_sum <- character(length(each_tumor_uniq_sig_variants))
#   gene_name_sum <- character(length(each_tumor_uniq_sig_variants))
#   #sum_target_samples_mut_number <- numeric(length(each_tumor_uniq_sig_variants))
#   #sum_target_samples_mut_callable <- numeric(length(each_tumor_uniq_sig_variants))
#   #sum_target_sample_mut_genome_mut_number <- numeric(length(each_tumor_uniq_sig_variants))
#   #sum_target_sample_mut_genome_mut_callable<- numeric(length(each_tumor_uniq_sig_variants))
#   
#   total_numbersamples_withgene <- 0
#   total_numbersamples_withoutgene <- 0 
#   
#   # calculate total x and total y
#   for (index in 1:length(each_tumor_uniq_sig_variants)) {
#     
#     print(paste("Processing the ", index, "variants, with variants", length(each_tumor_uniq_sig_variants),sep = " "))
#     each_variant <- each_tumor_uniq_sig_variants[index]
#     each_gene <- each_tumor_info[chrom_loc==each_variant]$gene_name[1]
#     each_ensembl <- each_tumor_info[chrom_loc==each_variant]$ensembl_id[1]
#     total_others <- 0
#     total_others_without <- 0
#     sample_loc <- which(each_variant == each_tumor_variants)
#     
#     numbersamples_withgene <- length(unique(each_tumor_info$sample_names[sample_loc]))
#     numbersamples_withoutgene <- each_tumor_total_sample- numbersamples_withgene
#     
#     total_numbersamples_withgene <- total_numbersamples_withgene + numbersamples_withgene
#     total_numbersamples_withoutgene <- total_numbersamples_withoutgene + numbersamples_withoutgene
#     
#     #target_samples_mut_number <-  sum(each_tumor_info[chrom_loc==each_variant,.(ensembl_mut_numer)])
#     #target_samples_mut_callable <-  sum(each_tumor_info[chrom_loc==each_variant,.(ensembl_callable)])
#     
#     #target_sample_mut_genome_mut_number <- sum(each_tumor_info[chrom_loc==each_variant,.(sample_genome_wide_mut_number)])
#     #target_sample_mut_genome_mut_callable <- sum(each_tumor_info[chrom_loc==each_variant,.(sample_genome_wide_mut_callable)])
#     ## fill up the info and create a table 
#     
#     gene_name_sum[index] <- each_gene
#     ensmbl_id_sum[index] <- each_ensembl
#     tumor_type_sum[index] <- each_tumor
#     gene_mut_sum[index] <- each_variant
#     numbersamples_withgene_sum[index] <- numbersamples_withgene
#     numbersamples_withoutgene_sum[index] <- numbersamples_withoutgene
#     
#     #sum_target_samples_mut_number[index] <- target_samples_mut_number
#     #sum_target_samples_mut_callable[index] <- target_samples_mut_callable
#     #sum_target_sample_mut_genome_mut_number[index] <- target_sample_mut_genome_mut_number
#     #sum_target_sample_mut_genome_mut_callable[index] <- target_sample_mut_genome_mut_callable
#     
#   }
#   total_others_sum <- total_numbersamples_withgene -numbersamples_withgene_sum
#   total_others_without_sum <- total_numbersamples_withoutgene- numbersamples_withoutgene_sum
#   
#   
#   
#   
#   each_tumor_info_sum <- data.table( tumor_type = tumor_type_sum,
#                                      ensembl_id =ensmbl_id_sum,
#                                      gene_name = gene_name_sum,
#                                      gene_mut = gene_mut_sum,
#                                      numbersamples_withgene = numbersamples_withgene_sum,
#                                      numbersamples_withoutgene = numbersamples_withoutgene_sum,
#                                      total_others = total_others_sum,
#                                      total_others_withoutgene = total_others_without_sum
#                                      #total_gene_mut_number = sum_target_samples_mut_number,
#                                      #total_gene_mut_callable = sum_target_samples_mut_callable,
#                                      #target_sample_genome_mut_number = sum_target_sample_mut_genome_mut_number,
#                                      #target_sample_genome_mut_callable = sum_target_sample_mut_genome_mut_callable
#                                      
#   )
#   
#   p_value_sum <- numeric(length(each_tumor_uniq_sig_variants))
#   for ( i in 1:nrow(each_tumor_info_sum)) {
#     
#     target <- as.matrix(each_tumor_info_sum[i,.(numbersamples_withgene,numbersamples_withoutgene)])
#     others <- as.matrix(each_tumor_info_sum[i,.(total_others,total_others_withoutgene)])
#     testor <- rbind(c(target),c(others))
#     each_gene_p_value <- fisher.test(testor, alternative = "greater")$p.value
#     p_value_sum[i] <-  each_gene_p_value
#     
#   }
#   each_tumor_info_sum$p_value <- p_value_sum 
#   each_tumor_info_sum <- each_tumor_info_sum[order(p_value)]
#   each_tumor_info_sum$BH_pvalue = p.adjust(each_tumor_info_sum$p_value, method = "BH")
#   
#   total_variant_summary <- rbindlist(list(total_variant_summary, each_tumor_info_sum))
# }
# fwrite(total_variant_summary,
#        file = paste(output_dir,"final_breed_variants_nonsyn_tumorwide_VAF_Mutect_orientBiasModified_04_07.txt",sep = seperator)
#        ,col.names = T,row.names = F,
#        quote = F,
#        eol = "\n",
#        sep ="\t")
# 


### Tumorwide genes
## assume all indel significant
sig_gene_sample_wide <- fread(paste(output_dir,"final_breeds_gene_nonsym_samplewide_p_value_orient_modify_04_07.txt",sep = seperator))
sig_gene_sample_wide <- sig_gene_sample_wide[,.(sample_names,gene_name,breeds_info,ensembl_id,Subtype,
                                                BH_pvalue)]

indel_file <- fread(paste(base_dir,"total_CDS_indel_info_withGene_04_08.txt",sep =seperator))
colnames(indel_file) <- c('chrom','pos','ref','alt','gene_name','ensembl_id','status','sample_names')
indel_file <- indel_file[gene_name!="-" & status!="nonframeshift " & ! sample_names %in% exclude,]
indel_file$Subtype <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table)
indel_file$breeds_info <- match_vector_table(indel_file$sample_names,"final_breed_label",whole_wes_clean_breed_table)

sig_indel_gene <- indel_file[,.(sample_names,gene_name,breeds_info,ensembl_id,Subtype)]
sig_indel_gene$BH_pvalue <- 0.04
sig_gene_sample_wide <- rbind(sig_gene_sample_wide,sig_indel_gene) 
sig_gene_sample_wide <- sig_gene_sample_wide[gene_name!="-",]

# fill up with NA string 
mutect_after_vaf <- fread(paste(base_dir,"Final_Total_withGene_Burair_Filtering3_VAF_Mutect_orientBiasModified_04_02.txt.gz",sep =seperator),
                          na.strings = "")
mutect_after_vaf <- mutect_after_vaf[status!= "synonymous",]
mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude]
mutect_after_vaf <- mutect_after_vaf[tumor_type!="UCL"]
mutect_after_vaf$Subtype <- match_vector_table(mutect_after_vaf$sample_names,column = "DiseaseAcronym2", table =whole_wes_clean_breed_table,string_value = T )
finalbreed <- match_vector_table(mutect_after_vaf$sample_names,column="final_breed_label",table=whole_wes_clean_breed_table)
mutect_after_vaf$Breeds <- finalbreed
mutect_after_vaf <- mutect_after_vaf[,chrom_loc:= paste(chrom,pos,sep = "_"),]

tumor_type <- unique(mutect_after_vaf$Subtype)
final_mutect_after_vaf <- mutect_after_vaf[,.(sample_names,gene_name,ensembl_id,status,ref,alt,Subtype)]
final_indel_file <- indel_file[,.(sample_names,gene_name,ensembl_id,status,ref,alt,Subtype)]
final_gene_tumor_wide_SNV_indel <- rbind(final_mutect_after_vaf,final_indel_file)
final_gene_tumor_wide_SNV_indel <- final_gene_tumor_wide_SNV_indel[!sample_names %in% exclude]

# 
# fwrite(final_gene_tumor_wide_SNV_indel,
#        file = paste(output_dir,"Combine_SNV_indel_Mutect_orientBiasModified_04_07.txt",sep = seperator)
#        ,col.names = T,row.names = F,
#        quote = F,
#        eol = "\n",
#        sep ="\t")


unique(final_gene_tumor_wide_SNV_indel[gene_name=="PTAFR" & Subtype=="TCL",.(gene_name,sample_names,status,ensembl_id)])
nrow(unique(final_gene_tumor_wide_SNV_indel[Subtype=="HSA",.(sample_names)]))
as.data.frame(table(sig_indel_gene[Subtype=="OM",.(gene_name)]))

### Tumorwide genes
total_gene_summary <- NULL
for (each_tumor in tumor_type) {
  current <- match(each_tumor, tumor_type,nomatch = 0)
  print(paste("Processing the ", current, "tumors with total numbers of tumor", length(tumor_type),sep = " "))
  each_tumor_info_sum <- NULL
  each_tumor_info <- final_gene_tumor_wide_SNV_indel[Subtype == each_tumor, ]
  each_tumor_genes <- each_tumor_info$gene_name
  each_tumor_sign_uniq_genes <- unique(sig_gene_sample_wide[BH_pvalue < 0.05 & Subtype==each_tumor, .(gene_name)]$gene_name)
  
  tumor_type_sum <- character(length(each_tumor_sign_uniq_genes))
  numbersamples_withgene_sum <-numeric(length(each_tumor_sign_uniq_genes))
  numbersamples_withoutgene_sum <- numeric(length(each_tumor_sign_uniq_genes))
  total_others_sum <- numeric(length(each_tumor_sign_uniq_genes))
  total_others_without_sum <- numeric(length(each_tumor_sign_uniq_genes))
  each_tumor_total_sample <-   length(unique(each_tumor_info$sample_names))
  ensmbl_id_sum <- character(length(each_tumor_sign_uniq_genes))
  gene_name_sum <- character(length(each_tumor_sign_uniq_genes))
  #sum_target_samples_mut_number <- numeric(length(each_tumor_sign_uniq_genes))
  #sum_target_samples_mut_callable <- numeric(length(each_tumor_sign_uniq_genes))
  #sum_target_sample_mut_genome_mut_number <- numeric(length(each_tumor_sign_uniq_genes))
  #sum_target_sample_mut_genome_mut_callable<- numeric(length(each_tumor_sign_uniq_genes))
  
  total_numbersamples_withgene <- 0
  total_numbersamples_withoutgene <- 0 
  
  # calculate total x and total y
  for (index in 1:length(each_tumor_sign_uniq_genes)) {
    #current <- match(each_gene, each_tumor_genes,nomatch = 0)
    print(paste("Processing the ", index, "genes, with total_genes", length(each_tumor_sign_uniq_genes),sep = " "))
    each_gene <- each_tumor_sign_uniq_genes[index]
    each_ensembl <- each_tumor_info[gene_name==each_gene]$ensembl_id[1]
    total_others <- 0
    total_others_without <- 0
    sample_loc <- which(each_tumor_genes==each_gene)
    
    #length(unique(each_tumor_info[gene_mutation == each_gene, ]$sample_names))
    numbersamples_withgene <- length(unique(each_tumor_info$sample_names[sample_loc]))
    numbersamples_withoutgene <- each_tumor_total_sample- numbersamples_withgene
    
    total_numbersamples_withgene <- total_numbersamples_withgene+numbersamples_withgene
    total_numbersamples_withoutgene <- total_numbersamples_withoutgene+numbersamples_withoutgene
    
    #target_samples_mut_number <-  sum(each_tumor_info[gene_name==each_gene,.(ensembl_mut_numer)])
    #target_samples_mut_callable <-  sum(each_tumor_info[gene_name==each_gene,.(ensembl_callable)])
    
    #target_sample_mut_genome_mut_number <- sum(each_tumor_info[gene_name==each_gene,.(sample_genome_wide_mut_number)])
    #target_sample_mut_genome_mut_callable <- sum(each_tumor_info[gene_name==each_gene,.(sample_genome_wide_mut_callable)])
    
    gene_name_sum[index] <- each_gene
    ensmbl_id_sum[index] <- each_ensembl
    tumor_type_sum[index] <- each_tumor
    numbersamples_withgene_sum[index] <- numbersamples_withgene
    numbersamples_withoutgene_sum[index] <- numbersamples_withoutgene
    #sum_target_samples_mut_number[index] <- target_samples_mut_number
    #sum_target_samples_mut_callable[index] <- target_samples_mut_callable
    #sum_target_sample_mut_genome_mut_number[index] <- target_sample_mut_genome_mut_number
    #sum_target_sample_mut_genome_mut_callable[index] <- target_sample_mut_genome_mut_callable
    
  }
  total_others_sum <- total_numbersamples_withgene -numbersamples_withgene_sum
  total_others_without_sum <- total_numbersamples_withoutgene- numbersamples_withoutgene_sum
  
  
  
  
  each_tumor_info_sum <- data.table( tumor_type = tumor_type_sum,
                                     ensembl_id =ensmbl_id_sum,
                                     gene_name = gene_name_sum,
                                     numbersamples_withgene = numbersamples_withgene_sum,
                                     numbersamples_withoutgene = numbersamples_withoutgene_sum,
                                     total_others = total_others_sum,
                                     total_others_withoutgene = total_others_without_sum
                                     #total_gene_mut_number = sum_target_samples_mut_number,
                                     #total_gene_mut_callable = sum_target_samples_mut_callable,
                                     #target_sample_genome_mut_number = sum_target_sample_mut_genome_mut_number,
                                     #target_sample_genome_mut_callable = sum_target_sample_mut_genome_mut_callable)
  )
  
  p_value_sum <- numeric(length(each_tumor_sign_uniq_genes))
  for ( i in 1:nrow(each_tumor_info_sum)) {
    
    target <- as.matrix(each_tumor_info_sum[i,.(numbersamples_withgene,numbersamples_withoutgene)])
    others <- as.matrix(each_tumor_info_sum[i,.(total_others,total_others_withoutgene)])
    testor <- rbind(c(target),c(others))
    each_gene_p_value <- fisher.test(testor, alternative = "greater")$p.value
    p_value_sum[i] <-  each_gene_p_value
    
  }
  each_tumor_info_sum$p_value <- p_value_sum 
  each_tumor_info_sum <- each_tumor_info_sum[order(p_value)]
  each_tumor_info_sum$BH_pvalue = p.adjust(each_tumor_info_sum$p_value, method = "BH")
  
  total_gene_summary <- rbindlist(list(total_gene_summary, each_tumor_info_sum))
}

top_8_tumor_gene <- total_gene_summary[, head(.SD,8 ), by=.(tumor_type)]

fwrite(top_8_tumor_gene,
       file = paste(output_dir,"Top8_tumorwide_VAF_Mutect_orientBiasModified_04_09.txt",sep = seperator)
       ,col.names = T,row.names = F,
       quote = F,
       eol = "\n",
       sep ="\t")


fwrite(total_gene_summary,
       file = paste(output_dir,"final_breed_genes_nonsyn_tumorwide_VAF_Mutect_orientBiasModified_04_07.txt",sep = seperator)
       ,col.names = T,row.names = F,
       quote = F,
       eol = "\n",
       sep ="\t")

