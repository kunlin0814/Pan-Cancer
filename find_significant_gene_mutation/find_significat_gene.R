library(data.table)
library(tidyverse)
library(readxl)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
  #"/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/Burair_filtering3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1/"
seperator <- "/"

#retro_gene_list <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Retro_gene_finding/RetroGeneList/new_retro_gene_list_CanFam3.1.99gtf.txt",
#                         header = F)
whole_wes_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt") 
             #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt")       

exclude <- unique(unlist(whole_wes_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))


## check duplicated
# a <- gene[sample_names=="CCB040105"& ensembl_id=="ENSCAFG00000001781"]

# 
Breed_info <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan-Cancer-Manuscript/Figure1/WES_WGS_merge.xlsx",
                         sheet ="WES_WGS")

Breed_info <- setDT(Breed_info)
mutect_after_vaf <- fread(paste(base_dir,"02_01","Burair_WithBreeds_Subtypes_QCpass_filtering3_mutect_after_vaf_02_01.txt",sep =seperator))

# subtype <- sapply(table_total_sample,FUN = match_table, column="DiseaseAcronym2",table=whole_wes_table)
# breed <- sapply(table_total_sample,FUN = match_table, column="Breeds",table=Breed_info)
# mutect_after_vaf$Subtype <- subtype
# mutect_after_vaf$Breeds <- breed
# mutect_after_vaf <- mutect_after_vaf[,chrom_loc:= paste(chrom,pos,sep = "_"),]
# mutect_after_vaf$gene_TMB <- (mutect_after_vaf$ensembl_mut_numer*1000000)/mutect_after_vaf$ensembl_callable
# mutect_after_vaf$genome_TMB <- (mutect_after_vaf$sample_genome_wide_mut_number*1000000)/mutect_after_vaf$sample_genome_wide_mut_callable
# mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude & gene_name!="-", ]

## get_rid of syn mut
mutect_after_vaf <- mutect_after_vaf[status!= "synonymous",]

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
       file = paste(base_dir,"02_01","variant_nonsyn_samplewide_p_value_VAF_Mutect_orientBias3_02_01.gz",sep = seperator)
       ,col.names = T,
       row.names = F,
       quote = F,
       eol = "\n",
       compress = "gzip",
       sep ="\t")


### samplewide gene_name ##

gene_total_info_sum <- NULL
for (index in 1:length(total_sample)) {
  print(paste("processing the",index,"sample, with total samples",length(total_sample),sep = " " ))
  sample <- total_sample[index]
  info_sum <- NULL
  subtype <-  mutect_after_vaf[sample_names==sample,.(Subtype)]$Subtype[1]
  tumor <-  mutect_after_vaf[sample_names==sample,.(tumor_type)]$tumor_type[1]
  gene_name <- mutect_after_vaf[sample_names==sample,.(gene_name)]
  each_sample_gene_id <- sort(unique(gene_name[["gene_name"]]))
  for (i in each_sample_gene_id ){
    ensembl_id <- mutect_after_vaf[sample_names==sample & gene_name==i, .(ensembl_id)]$ensembl_id[1]
    target <- mutect_after_vaf[sample_names==sample & gene_name==i, .(tRef,tAlt)]
    target_combine <- target[, .(tRef = sum(tRef), tAlt = sum(tAlt)),]
    
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
                       ensembl_id = ensembl_id,
                       tumor_type = tumor,
                       subtype = subtype,
                       p_value = p_value,
                       ensembl_mut_number =ensembl_mut_number,
                       ensembl_callable = ensembl_callable,
                       sample_genome_wide_mut_number = sample_genome_wide_mut_number,
                       sample_genome_wide_mut_callable = sample_genome_wide_mut_callable
                       )
    info_sum <- rbindlist(list(info_sum,info))
  }
  info_sum <- info_sum[order(p_value)]
  info_sum$BH_pvalue = p.adjust(info_sum$p_value, method = "BH")
  gene_total_info_sum <- rbindlist(list(gene_total_info_sum,info_sum))
}

gene_total_info_sum$gene_TMB <- (gene_total_info_sum$ensembl_mut_number*1000000)/gene_total_info_sum$ensembl_callable
gene_total_info_sum$genome_TMB <- (gene_total_info_sum$sample_genome_wide_mut_number*1000000)/gene_total_info_sum$sample_genome_wide_mut_callable



fwrite(gene_total_info_sum,
       file = paste(base_dir,"02_01","gene_nonsym_samplewide_p_value_VAF_Mutect_orientBias3_02_01.gz",sep = seperator)
       ,col.names = T,row.names = F,
       quote = F,
       eol = "\n",
       compress = "gzip",
       sep ="\t")


#a <- which("chr34_12675674" == 
MT <- mutect_after_vaf[tumor_type=="MT",]

#mutect_after_vaf$sample_names[sample_loc])

sig_variants_sample_wide <- fread(file = paste(base_dir,"02_01","variant_nonsyn_samplewide_p_value_VAF_Mutect_orientBias3_02_01.gz",sep = seperator))

sig_gene_sample_wide <- fread(paste(base_dir,"02_01","gene_nonsym_samplewide_p_value_VAF_Mutect_orientBias3_02_01.gz",sep = seperator))


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
  
    target_samples_mut_number <-  sum(each_tumor_info[chrom_loc==each_variant,.(ensembl_mut_number)])
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
                                     total_others_withoutgene = total_others_without_sum,
                                     total_gene_mut_number = sum_target_samples_mut_number,
                                     total_gene_mut_callable = sum_target_samples_mut_callable,
                                     target_sample_genome_mut_number = sum_target_sample_mut_genome_mut_number,
                                     target_sample_genome_mut_callable = sum_target_sample_mut_genome_mut_callable
                                     
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
       file = paste(base_dir,"02_01","variant_nonsyn_tumorwide_p_value_VAF_Mutect_orientBias3_02_01.gz",sep = seperator)
       ,col.names = T,row.names = F,
       eol = "\n",
       quote = F,
       compress = "gzip",
       sep ="\t")


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
    
    target_samples_mut_number <-  sum(each_tumor_info[chrom_loc==each_gene,.(ensembl_mut_number)])
    target_samples_mut_callable <-  sum(each_tumor_info[chrom_loc==each_gene,.(ensembl_callable)])
    
    target_sample_mut_genome_mut_number <- sum(each_tumor_info[chrom_loc==each_gene,.(sample_genome_wide_mut_number)])
    target_sample_mut_genome_mut_callable <- sum(each_tumor_info[chrom_loc==each_gene,.(sample_genome_wide_mut_callable)])
    
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
                                     total_others_withoutgene = total_others_without_sum,
                                     total_gene_mut_number = sum_target_samples_mut_number,
                                     total_gene_mut_callable = sum_target_samples_mut_callable,
                                     target_sample_genome_mut_number = sum_target_sample_mut_genome_mut_number,
                                     target_sample_genome_mut_callable = sum_target_sample_mut_genome_mut_callable)
  
  
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
       file = paste(base_dir,"02_01","gene_nonsyn_tumorwide_p_value_VAF_Mutect_orientBias3_02_01.gz",sep = seperator)
      ,col.names = T,row.names = F,
       quote = F,
      eol = "\n",
      compress = "gzip",
      sep ="\t")

# 
# ### check the results

source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")
base_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/Burair_filtering3/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/VAF/New_Burair_filterin3/Mutect1/"
seperator <- "/"


variant_sample <- fread(paste(base_dir,"02_01","variant_nonsyn_samplewide_p_value_VAF_Mutect_orientBias3_02_01.gz",
                              sep = seperator))
variant_tumor <- fread(paste(base_dir,"02_01","variant_nonsyn_tumorwide_p_value_VAF_Mutect_orientBias3_02_01.gz",
                             sep = seperator))
gene_sample <- fread(paste(base_dir,"02_01","gene_nonsyn_samplewide_p_value_VAF_Mutect_orientBias3_02_01.gz",
                           sep = seperator))
gene_tumor <- fread(paste(base_dir,"02_01","gene_nonsyn_tumorwide_p_value_VAF_Mutect_orientBias3_02_01.gz",
                          sep = seperator))

mutect_after_vaf[gene_name=="NRAS" & tumor_type=="HSA",.(sample_names)]

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
