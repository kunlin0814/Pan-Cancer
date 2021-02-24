library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
  #"/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- 
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1/"
seperator <- "/"

#retro_gene_list <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Retro_gene_finding/RetroGeneList/new_retro_gene_list_CanFam3.1.99gtf.txt",
#                         header = F)
## check duplicated
# a <- gene[sample_names=="CCB040105"& ensembl_id=="ENSCAFG00000001781"]

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 



exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))


# fill up with NA string 
mutect_after_vaf <- fread(paste(base_dir,"02_11","mutect_noucl_vaf_withBreeds_callable_0210.txt",sep =seperator),
                          na.strings = "")

mutect_after_vaf <- mutect_after_vaf[status!= "synonymous",]
mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude]

 

# fwrite(mutect_after_vaf, file = paste(base_dir,"02_18","NonSyn_Burair_filtering3_WithBreeds_Subtypes_QCpass_mutect_after_vaf_02_18.txt",sep = seperator),
#        col.names = T, row.names = F, quote = F, sep = "\t")


### append column ###
# subtype <- match_vector_table(mutect_after_vaf$sample_names,column = "DiseaseAcronym2", table =whole_wes_clean_breed_table,string_value = T )
# mutect_after_vaf$Subtype <- subtype

finalbreed <- match_vector_table(mutect_after_vaf$sample_names,column="final_breed_label",table=whole_wes_clean_breed_table)
mutect_after_vaf$Breeds <- finalbreed

# mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude & gene_name!="-", ]
# mutect_after_vaf <- mutect_after_vaf[,chrom_loc:= paste(chrom,pos,sep = "_"),]
# 

### append column end ###

total_sample <- unique( mutect_after_vaf$sample_names)
### samplewide variants ##
total_info_sum <- NULL
for (index in 1:length(total_sample)) {
  print(paste("processing the",index,"sample, with total samples",length(total_sample),sep = " " ))
  sample <- total_sample[index]
  info_sum <- NULL
  variant_loc <- mutect_after_vaf[sample_names==sample,.(chrom_loc)]
  for (i in unique(variant_loc[["chrom_loc"]]) ){
    info <- mutect_after_vaf[sample_names==sample & chrom_loc==i,]
    target <- mutect_after_vaf[sample_names==sample & chrom_loc==i, .(tRef,tAlt)]
    target_combine <- target[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    # other samples chrom_loc!=i
    others <- mutect_after_vaf[sample_names==sample & chrom_loc!=i, .(tRef,tAlt)]
    other_combine <- others[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    
    testor=rbindlist(list(target_combine,other_combine))
    p_value <- fisher.test(testor,alternative = "less")$p.value
    info <- info[,p_value:=p_value]
    info_sum <- rbindlist(list(info_sum,info))
  }
  info_sum <- info_sum[order(p_value)]
  info_sum$BH_pvalue = p.adjust(info_sum$p_value, method = "BH")
  total_info_sum <- rbindlist(list(total_info_sum,info_sum))
}

fwrite(total_info_sum,
       file = paste(base_dir,"02_18","final_breed_variant_nonsyn_samplewide_p_value_VAF_Mutect_orientBias3_02_18.txt",sep = seperator)
       ,col.names = T,
       row.names = F,
       quote = F,
       na = "NA",
       eol = "\n",
       sep ="\t")


### samplewide gene_name ##
total_sample <- unique( mutect_after_vaf$sample_names)
gene_total_info_sum <- NULL
for (index in 1:length(total_sample)) {
  print(paste("processing the",index,"sample, with total samples",length(total_sample),sep = " " ))
  sample <- total_sample[index]
  info_sum <- NULL
  subtype <-  mutect_after_vaf[sample_names==sample,.(Subtype)]$Subtype[1]
  tumor <-  mutect_after_vaf[sample_names==sample,.(tumor_type)]$tumor_type[1]
  gene_name <- mutect_after_vaf[sample_names==sample,.(gene_name)]
  each_sample_gene_id <- sort(unique(gene_name[["gene_name"]]))
  breed <- mutect_after_vaf[sample_names==sample,.(Breeds)]$Breeds[1]
  for (i in each_sample_gene_id ){
    ensembl_id <- mutect_after_vaf[sample_names==sample & gene_name==i, .(ensembl_id)]$ensembl_id[1]
    target <- mutect_after_vaf[sample_names==sample & gene_name==i, .(tRef,tAlt)]
    target_combine <- target[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    mut_status <- mutect_after_vaf[sample_names==sample & gene_name==i, .(status)]$status[1]
    # others gene_name!=i
    others <- mutect_after_vaf[sample_names==sample & gene_name!=i, .(tRef,tAlt)]
    other_combine <- others[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    
    testor=rbindlist(list(target_combine,other_combine))
    p_value <- fisher.test(testor,alternative = "less")$p.value
    
    ensembl_callable <- mutect_after_vaf[sample_names==sample & ensembl_id==ensembl_id, .(ensembl_callable)]$ensembl_callable[1]
    ensembl_mut_number <- mutect_after_vaf[sample_names==sample & ensembl_id==ensembl_id, .(ensembl_mut_numer)]$ensembl_mut_numer[1]
    sample_genome_wide_mut_number <- mutect_after_vaf[sample_names==sample & ensembl_id==ensembl_id, .(sample_genome_wide_mut_number)]$sample_genome_wide_mut_number[1]
    sample_genome_wide_mut_callable <- mutect_after_vaf[sample_names==sample & ensembl_id==ensembl_id, .(sample_genome_wide_mut_callable)]$sample_genome_wide_mut_callable[1]
    
    info <- data.table(sample_names = sample, 
                       gene_name= i, 
                       breeds_info = breed,
                       ensembl_id = ensembl_id,
                       tumor_type = tumor,
                       subtype = subtype,
                       p_value = p_value,
                       #ensembl_mut_number =ensembl_mut_number,
                       #ensembl_callable = ensembl_callable,
                       #sample_genome_wide_mut_number = sample_genome_wide_mut_number,
                       #sample_genome_wide_mut_callable = sample_genome_wide_mut_callable,
                       mut_type = mut_status
                       )
    info_sum <- rbindlist(list(info_sum,info))
  }
  info_sum <- info_sum[order(p_value)]
  info_sum$BH_pvalue = p.adjust(info_sum$p_value, method = "BH")
  gene_total_info_sum <- rbindlist(list(gene_total_info_sum,info_sum))
}

#gene_total_info_sum$gene_TMB <- (gene_total_info_sum$ensembl_mut_number*1000000)/gene_total_info_sum$ensembl_callable
#gene_total_info_sum$genome_TMB <- (gene_total_info_sum$sample_genome_wide_mut_number*1000000)/gene_total_info_sum$sample_genome_wide_mut_callable

# taifang_data <- gene_total_info_sum[,.(sample_names,gene_name,tumor_type,mut_type)]


fwrite(gene_total_info_sum,
       file = paste(base_dir,"02_18","final_breeds_gene_nonsym_samplewide_p_value_VAF_Mutect_orientBias3_02_18.txt",sep = seperator)
       ,col.names = T,row.names = F,
       quote = F,
       eol = "\n",
       sep ="\t",
       na = "NA")

#### Tumor wide analysis ####

sig_variants_sample_wide <- fread(file = paste(base_dir,"02_18","final_breed_variant_nonsyn_samplewide_p_value_VAF_Mutect_orientBias3_02_18.txt",sep = seperator))

sig_gene_sample_wide <- fread(paste(base_dir,"02_18","final_breeds_gene_nonsym_samplewide_p_value_VAF_Mutect_orientBias3_02_18.txt",sep = seperator))


whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))




# fill up with NA string 
mutect_after_vaf <- fread(paste(base_dir,"02_11","mutect_noucl_vaf_withBreeds_callable_0210.txt",sep =seperator),
                          na.strings = "")

mutect_after_vaf <- mutect_after_vaf[status!= "synonymous",]
mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude]
## change breed info into final breed info
finalbreed <- match_vector_table(mutect_after_vaf$sample_names,column="final_breed_label",table=whole_wes_clean_breed_table)
mutect_after_vaf$Breeds <- finalbreed

## Tumorwide
tumor_type <- unique(mutect_after_vaf$Subtype)

### Tumorwide variants

total_variant_summary <- NULL
for (each_tumor in tumor_type) {
  current <- match(each_tumor, tumor_type,nomatch = 0)
  print(paste("Processing the ", current, "tumors with total numbers of tumor", length(tumor_type),sep = " "))
  each_tumor_info_sum <- NULL
  each_tumor_info <- mutect_after_vaf[Subtype == each_tumor, ]
  each_tumor_variants <- each_tumor_info$chrom_loc
  
  # the sig variants come from samplewide vaf comparision 
  each_tumor_uniq_sig_variants <- unique(sig_variants_sample_wide[BH_pvalue < 0.05 & Subtype==each_tumor, .(chrom_loc)]$chrom_loc)
  
  ## vectorized the value I want to append later
  tumor_type_sum <- character(length(each_tumor_uniq_sig_variants))
  gene_mut_sum <- character(length(each_tumor_uniq_sig_variants))
  numbersamples_withgene_sum <-numeric(length(each_tumor_uniq_sig_variants))
  numbersamples_withoutgene_sum <- numeric(length(each_tumor_uniq_sig_variants))
  total_others_sum <- numeric(length(each_tumor_uniq_sig_variants))
  total_others_without_sum <- numeric(length(each_tumor_uniq_sig_variants))
  each_tumor_total_sample <-   length(unique(each_tumor_info$sample_names))
  ensmbl_id_sum <- character(length(each_tumor_uniq_sig_variants))
  gene_name_sum <- character(length(each_tumor_uniq_sig_variants))
  sum_target_samples_mut_number <- numeric(length(each_tumor_uniq_sig_variants))
  sum_target_samples_mut_callable <- numeric(length(each_tumor_uniq_sig_variants))
  sum_target_sample_mut_genome_mut_number <- numeric(length(each_tumor_uniq_sig_variants))
  sum_target_sample_mut_genome_mut_callable<- numeric(length(each_tumor_uniq_sig_variants))
  
  total_numbersamples_withgene <- 0
  total_numbersamples_withoutgene <- 0 
  
  # calculate total x and total y
  for (index in 1:length(each_tumor_uniq_sig_variants)) {
    
    print(paste("Processing the ", index, "variants, with variants", length(each_tumor_uniq_sig_variants),sep = " "))
    each_variant <- each_tumor_uniq_sig_variants[index]
    each_gene <- each_tumor_info[chrom_loc==each_variant]$gene_name[1]
    each_ensembl <- each_tumor_info[chrom_loc==each_variant]$ensembl_id[1]
    total_others <- 0
    total_others_without <- 0
    sample_loc <- which(each_variant == each_tumor_variants)
    
    numbersamples_withgene <- length(unique(each_tumor_info$sample_names[sample_loc]))
    numbersamples_withoutgene <- each_tumor_total_sample- numbersamples_withgene
    
    total_numbersamples_withgene <- total_numbersamples_withgene + numbersamples_withgene
    total_numbersamples_withoutgene <- total_numbersamples_withoutgene + numbersamples_withoutgene
  
    target_samples_mut_number <-  sum(each_tumor_info[chrom_loc==each_variant,.(ensembl_mut_numer)])
    target_samples_mut_callable <-  sum(each_tumor_info[chrom_loc==each_variant,.(ensembl_callable)])
    
    target_sample_mut_genome_mut_number <- sum(each_tumor_info[chrom_loc==each_variant,.(sample_genome_wide_mut_number)])
    target_sample_mut_genome_mut_callable <- sum(each_tumor_info[chrom_loc==each_variant,.(sample_genome_wide_mut_callable)])
    ## fill up the info and create a table 
    
    gene_name_sum[index] <- each_gene
    ensmbl_id_sum[index] <- each_ensembl
    tumor_type_sum[index] <- each_tumor
    gene_mut_sum[index] <- each_variant
    numbersamples_withgene_sum[index] <- numbersamples_withgene
    numbersamples_withoutgene_sum[index] <- numbersamples_withoutgene

    sum_target_samples_mut_number[index] <- target_samples_mut_number
    sum_target_samples_mut_callable[index] <- target_samples_mut_callable
    sum_target_sample_mut_genome_mut_number[index] <- target_sample_mut_genome_mut_number
    sum_target_sample_mut_genome_mut_callable[index] <- target_sample_mut_genome_mut_callable
    
  }
  total_others_sum <- total_numbersamples_withgene -numbersamples_withgene_sum
  total_others_without_sum <- total_numbersamples_withoutgene- numbersamples_withoutgene_sum
  
  
  
  
  each_tumor_info_sum <- data.table( tumor_type = tumor_type_sum,
                                     ensembl_id =ensmbl_id_sum,
                                     gene_name = gene_name_sum,
                                     gene_mut = gene_mut_sum,
                                     numbersamples_withgene = numbersamples_withgene_sum,
                                     numbersamples_withoutgene = numbersamples_withoutgene_sum,
                                     total_others = total_others_sum,
                                     total_others_withoutgene = total_others_without_sum
                                     #total_gene_mut_number = sum_target_samples_mut_number,
                                     #total_gene_mut_callable = sum_target_samples_mut_callable,
                                     #target_sample_genome_mut_number = sum_target_sample_mut_genome_mut_number,
                                     #target_sample_genome_mut_callable = sum_target_sample_mut_genome_mut_callable
                                     
  )
  
  p_value_sum <- numeric(length(each_tumor_uniq_sig_variants))
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
  
  total_variant_summary <- rbindlist(list(total_variant_summary, each_tumor_info_sum))
}


fwrite(total_variant_summary,
       file = paste(base_dir,"02_18","final_breed_variant_nonsyn_tumorwide_p_value_VAF_Mutect_orientBias3_02_18.txt",sep = seperator)
       ,col.names = T,row.names = F,
       eol = "\n",
       quote = F,
       sep ="\t",
       na = "NA")


## Tumorwide
## 2/2
## use mutect_vaf input file to identify the occurence 
## for variants and genes, if the variants or gene are not sig if all samples, then the variants and gene won't be choosed
## The variants or genes must be sig in at least one samples, so that it will be choosed 


tumor_type <- unique(mutect_after_vaf$Subtype)
### Tumorwide genes
total_gene_summary <- NULL
for (each_tumor in tumor_type) {
  current <- match(each_tumor, tumor_type,nomatch = 0)
  print(paste("Processing the ", current, "tumors with total numbers of tumor", length(tumor_type),sep = " "))
  each_tumor_info_sum <- NULL
  each_tumor_info <- mutect_after_vaf[Subtype == each_tumor, ]
  each_tumor_genes <- each_tumor_info$gene_name
  each_tumor_sign_uniq_genes <- unique(sig_gene_sample_wide[BH_pvalue < 0.05 & subtype==each_tumor, .(gene_name)]$gene_name)
  
  tumor_type_sum <- character(length(each_tumor_sign_uniq_genes))
  numbersamples_withgene_sum <-numeric(length(each_tumor_sign_uniq_genes))
  numbersamples_withoutgene_sum <- numeric(length(each_tumor_sign_uniq_genes))
  total_others_sum <- numeric(length(each_tumor_sign_uniq_genes))
  total_others_without_sum <- numeric(length(each_tumor_sign_uniq_genes))
  each_tumor_total_sample <-   length(unique(each_tumor_info$sample_names))
  ensmbl_id_sum <- character(length(each_tumor_sign_uniq_genes))
  gene_name_sum <- character(length(each_tumor_sign_uniq_genes))
  sum_target_samples_mut_number <- numeric(length(each_tumor_sign_uniq_genes))
  sum_target_samples_mut_callable <- numeric(length(each_tumor_sign_uniq_genes))
  sum_target_sample_mut_genome_mut_number <- numeric(length(each_tumor_sign_uniq_genes))
  sum_target_sample_mut_genome_mut_callable<- numeric(length(each_tumor_sign_uniq_genes))
  
  
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
    sample_loc <- which(each_gene == each_tumor_genes)
    
    #length(unique(each_tumor_info[gene_mutation == each_gene, ]$sample_names))
    numbersamples_withgene <- length(unique(each_tumor_info$sample_names[sample_loc]))
    numbersamples_withoutgene <- each_tumor_total_sample- numbersamples_withgene
    
    total_numbersamples_withgene <- total_numbersamples_withgene+numbersamples_withgene
    total_numbersamples_withoutgene <- total_numbersamples_withoutgene+numbersamples_withoutgene
    
    target_samples_mut_number <-  sum(each_tumor_info[gene_name==each_gene,.(ensembl_mut_numer)])
    target_samples_mut_callable <-  sum(each_tumor_info[gene_name==each_gene,.(ensembl_callable)])
    
    target_sample_mut_genome_mut_number <- sum(each_tumor_info[gene_name==each_gene,.(sample_genome_wide_mut_number)])
    target_sample_mut_genome_mut_callable <- sum(each_tumor_info[gene_name==each_gene,.(sample_genome_wide_mut_callable)])
    
    gene_name_sum[index] <- each_gene
    ensmbl_id_sum[index] <- each_ensembl
    tumor_type_sum[index] <- each_tumor
    numbersamples_withgene_sum[index] <- numbersamples_withgene
    numbersamples_withoutgene_sum[index] <- numbersamples_withoutgene
    sum_target_samples_mut_number[index] <- target_samples_mut_number
    sum_target_samples_mut_callable[index] <- target_samples_mut_callable
    sum_target_sample_mut_genome_mut_number[index] <- target_sample_mut_genome_mut_number
    sum_target_sample_mut_genome_mut_callable[index] <- target_sample_mut_genome_mut_callable
    
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


fwrite(total_gene_summary,
       file = paste(base_dir,"02_18","final_breed_gene_nonsyn_tumorwide_p_value_VAF_Mutect_orientBias3_02_18.txt",sep = seperator)
      ,col.names = T,row.names = F,
       quote = F,
      eol = "\n",
      sep ="\t")

#### Breeds within each tumor type 
# need at least 10 dogs,
# need at least two certain dogs have that mutation (gene or variants)
base_dir <- 
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1/"
seperator <- "/"
whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))


SNV <- fread(paste(base_dir,"02_11","mutect_noucl_vaf_withBreeds_callable_0210.txt",sep =seperator)
                          ,fill = T,na.strings="")
SNV <- SNV[status!= "synonymous",]
SNV <- SNV[!sample_names %in% exclude,]
finalbreed <- match_vector_table(SNV$sample_names,column="final_breed_label",table=whole_wes_clean_breed_table)
SNV$Breeds <- finalbreed

# total_breeds <- unique( SNV$Breeds)
# clean_breeds <- na.omit(total_breeds)
indel_file <- fread(paste(base_dir,"passQC_pan-tumor-total_indel_info_0214.txt",sep =seperator))
indel_file <- indel_file[gene_name!="-" & status=="frameshift" & ! sample_names %in% exclude,]
setcolorder(indel_file,c("sample_names","gene_name","emsembl_id","status"))
indel_file <- indel_file[,emsembl_id:=NULL]

Subtype <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table )
indel_file$Subtype <- Subtype
breed <- match_vector_table(indel_file$sample_names,"final_breed_label",whole_wes_clean_breed_table )
indel_file$Breeds <- breed
final_SNV <- SNV[,.(sample_names,gene_name,status,Subtype,Breeds)]

mutect_after_vaf <- rbindlist(list(final_SNV,indel_file))

mutect_after_vaf <- mutect_after_vaf[Subtype!="UCL" & !sample_names %in% exclude,]

#### Breeds within each tumor type 
# need at least 10 dogs,
# need at least two certain dogs have that mutation (gene or variants)

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
                       !Breeds %in% c("Mixed",NA)]
  if (nrow(each_tumor_info) > 0) {
    ## check breeds number
    each_tumor_all_breeds <-
      unique(each_tumor_info[, .(sample_names, Breeds)])
    each_tumor_breed_number <-
      data.table(table(each_tumor_all_breeds$Breeds))
    candidate_breeds <-
      sort(each_tumor_breed_number[N >= number_breeds_cutoff,]$V1)
    #candidate_breeds <- candidate_breeds[!candidate_breeds %in% c("Mixed")]
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

meet_cut_off <- total_tumor_type_summary[target_breeds_with >number_sample_mut_cutoff,]

## within each tumor type, do p adjustment for each breed

Total_tumor_info <- NULL
for (each_tumor_type in tumor_type){
  each_tumor_final_breed_label <- NULL
  #each_tumor_type <- "MT"
  candidate_breed_each_tumor_type <- unique(unlist(meet_cut_off[tumor_type ==each_tumor_type,.(breeds)]$breeds))
  for ( each_breed in candidate_breed_each_tumor_type){
    each_breed_pvalue_for_each_tumor <- meet_cut_off[tumor_type == each_tumor_type & breeds == each_breed]
    each_breed_pvalue_for_each_tumor <- each_breed_pvalue_for_each_tumor[order(p_value)]
    each_breed_pvalue_for_each_tumor$BH_pvalue = p.adjust(each_breed_pvalue_for_each_tumor$p_value, method = "BH")
    each_tumor_final_breed_label <- rbindlist(list(each_tumor_final_breed_label, each_breed_pvalue_for_each_tumor))
  }
  Total_tumor_info <- rbindlist(list(Total_tumor_info,each_tumor_final_breed_label))
}

fwrite(Total_tumor_info, file = paste(base_dir,"02_18","final_breed_WithBH_breed_significant_Tumor_wide_02_18.txt",sep=seperator),
       col.names = T, row.names = F, quote = F, eol = "\n", na = "NA",
       sep = "\t")




# ### check the results

# source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1/"
seperator <- "/"


variant_sample <- fread(paste(base_dir,"02_18","variant_nonsyn_samplewide_p_value_VAF_Mutect_orientBias3_02_18.txt",
                              sep = seperator))
variant_tumor <- fread(paste(base_dir,"02_18","variant_nonsyn_tumorwide_p_value_VAF_Mutect_orientBias3_02_18.txt",
                             sep = seperator))
gene_sample <- fread(paste(base_dir,"02_18","gene_nonsyn_samplewide_p_value_VAF_Mutect_orientBias3_02_18.txt",
                           sep = seperator))
gene_tumor <- fread(paste(base_dir,"02_18","gene_nonsyn_tumorwide_p_value_VAF_Mutect_orientBias3_02_18.txt",
                          sep = seperator))



# top_10_sample_variants <- variant_sample[, head(.SD, 10), by=.(tumor_type)]
# top_10_sample_genes <- gene_sample[, head(.SD, 10), by=.(tumor_type)]
top_6_tumor_variants <- variant_tumor[, head(.SD, 5), by=tumor_type]
top_6_tumor_gene <- gene_tumor[, head(.SD, 5), by=tumor_type]

fwrite(top_6_tumor_variants, file = paste(base_dir,"02_18","top_6_tumor_variants.txt",sep=seperator),
       col.names = T, row.names = F, quote = F, eol = "\n",
       sep = "\t")


fwrite(top_6_tumor_gene, file = paste(base_dir,"02_18","top_6_tumor_genes.txt",sep=seperator),
       col.names = T, row.names = F, quote = F, eol = "\n",
       sep = "\t")


gene_tumor[tumor_type=="GLM" & gene_name=="PIK3CA"]

mutect_after_vaf[gene_name=="PIK3CA" & tumor_type=="GLM",.(sample_names)]

gene_sign <- gene_tumor[tumor_type=="HSA"]
pik3_sample_num <- variant_sample[gene_name=="PIK3CA" & p_value<0.05, .(p_value,BH_pvalue), keyby = .(tumor_type,sample_names)]
tp53_sample_num <- variant_sample[gene_name=="TP53" & p_value<0.05, .N, keyby = .(tumor_type)]



tp53_sign <- unique(variant_sample[tumor_type=="HSA" & gene_name=="TP53",.(sample_names)])

freq <- as.data.frame(table(tumor$chrom_loc))

# file <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/Burair_filtering3/Mutect1/01_26/variant_tumorwide_p_value_Filtering3_VAF_Mutect_orientBias3.gz") 


# # 
# sample_variant <- fread(paste(base_dir,"significant","clean_variant_samplewide_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))
# sample_gene <- fread(paste(base_dir,"significant","clean_gene_samplewide_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))
# tumor_variant <- fread(paste(base_dir,"significant","clean_variant_tumorwide_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))
# tumor_gene <- fread(paste(base_dir,"significant","clean_gene_tumorwide_p_value_total_final_Filtering3_VAF_Mutect_orientBias3.gz",sep = seperator))
# 
# a <- sample_variant[sample_names=="004",]
# 
# significant_sample_variant <- sample_variant[!ensembl_id %in% retro_gene_list
#                                              $V1 & BH_pvalue < 0.05 
#                                              & gene_name!="-",.(chrom_loc, BH_pvalue,vaf), 
#                                              keyby = .(sample_names,tumor_type)]
# 
# significant_sample_gene <- sample_variant[!ensembl_id %in% retro_gene_list$V1 
#                                           & BH_pvalue < 0.05
#                                           & gene_name!="-",.(gene_name, BH_pvalue), 
#                                           keyby = .(sample_names,tumor_type)]
# significant_tumor_variant <- tumor_variant[!ensembl_id %in% retro_gene_list$V1 
#                                            & BH_pvalue < 0.05
#                                            & gene_name!="-",.(chrom_loc, BH_pvalue),
#                                            keyby = .(tumor_type)]
# significant_tumor_gene <- tumor_gene[!ensembl_id %in% retro_gene_list$V1 
#                                      & BH_pvalue < 0.05
#                                      & gene_name!="-",
#                                      .(gene_name, BH_pvalue), keyby = .(tumor_type)]
# 
# 
# top_10_sample_variants <- significant_sample_variant[, head(.SD, 10), by=.(sample_names)]
# top_10_sample_genes <- significant_sample_gene[, head(.SD, 10), by=.(sample_names)]
# top_10_tumor_variants <- significant_tumor_variant[, head(.SD, 10), by=tumor_type]
# top_10_tumor_gene <- significant_tumor_gene[, head(.SD, 10), by=tumor_type]
# 
# 
# # a <- top_10_sample_variants[tumor_type=="HSA",]
# # 
# # 
