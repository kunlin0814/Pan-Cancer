source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")

seperator <- "/"
target_tumor_type="OM"
human_base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_Cancer_paper/Human/OM"

pathway <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis/PathwayGeneList-03_19.xlsx",
                      sheet = "Pathway_gene_list")
pathway <- setDT(pathway)
target <- c("RTK/RAS","Cell cycle","PI3K","p53","Wnt","Chromatin remodeler")
target_pathway <- pathway[,target, with= F]

path_col <- colnames(target_pathway)


target_pathway_gene <- NULL

for (each_pathway in colnames(target_pathway)){
  each_path_gene <- target_pathway[, each_pathway, with = F][[each_pathway]]
  each_path_clean_gene <- each_path_gene[!is.na(each_path_gene)]
  target_pathway_gene <- c(target_pathway_gene,each_path_clean_gene)
}

### Dog data ###

base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_03_30.txt") 
#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_03_30.txt")

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))


## mutation Data

## SNV
mutect_after_vaf <- fread(paste(base_dir,"NonSyn_Burair_filtering3_WithBreeds_Subtypes_QCpass_mutect_after_vaf_02_11.txt",
                                sep =seperator))

# mutect_after_vaf <- fread(paste(base_dir,"total_final_withGene_final_Filtering3_VAF_Mutect1_orientBias3_0129.gz",
#                                 sep =seperator))
# mutect_after_vaf <- mutect_after_vaf[!sample_names %in% exclude]
# Subtype <- match_vector_table(mutect_after_vaf$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table )
# mutect_after_vaf$Subtype <- Subtype


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
#total_mut <- rbindlist(list(SNV,indel_file,CNV))
total_mut <- rbindlist(list(SNV,indel_file))

#total_sample <- unique(total_mut$sample_names)
total_sample <- unique(total_mut[Subtype=="OM",]$sample_names)

total_mut <- total_mut[!sample_names %in% exclude,]

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


## select target tumor type
## total_sum is pathway, total_mut is SNV+CNV
dog_target_tumor_pathway <- total_sum[Subtype==target_tumor_type,]

pathway_sum <- NULL
dog_total_sample_number <- nrow(dog_target_tumor_pathway)
for (each_pathway in colnames(target_pathway)){
  each_pathway_number <- sum(dog_target_tumor_pathway[,each_pathway, with =F ][[each_pathway]])
  each_sum <- list(Sample = each_pathway,
                   Mut_type = "path",
                   tumorxsAfflicted = each_pathway_number,
                   tumorsTotal = dog_total_sample_number)
  
  pathway_sum <- rbindlist(list(pathway_sum,each_sum))
}
pathway_sum <- setDT(pathway_sum)


#########  analyzed Human Data #########  
# sanger only has snv
#1. find the most frequeny mut gene for snv and cna in human 
base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_Cancer_paper/Human/OM"
human_snv <- fread(paste(base_dir,"sanger_human_data.txt",sep = seperator))
human_total_sample_number <- length(unique(human_snv$Sample))


total_snv_gene <- unique(human_snv$SYMBOL)

## count SNV
all_sample_snv_sum <- list()
for (each_gene in total_snv_gene){
  each_sample_number <- nrow(unique(human_snv[SYMBOL == each_gene,.(Sample,SYMBOL)])) 
  ratio = each_sample_number/human_total_sample_number
  each_gene_sum <- list(gene_name = each_gene,
                        numberOfsample = each_sample_number,
                        ratio = ratio)
  all_sample_snv_sum <-rbindlist(list(all_sample_snv_sum, each_gene_sum))
  
}
all_sample_snv_sum <- setDT(all_sample_snv_sum)
# 
# # count snv end 
# ## prepare CNV sample data
# human_cnv <- read.table("clipboard",sep = "\t",header = T,stringsAsFactors = F)
# human_cnv <- setDT(human_cnv)
# all_sample_cna <- list()
# # recreate Data format for cnv
# for (i in 1: nrow(human_cnv)){
#   gene_list <- unlist(strsplit(human_cnv[i,][['Gene']],","))
#   cn_status <- ifelse(human_cnv[i,][["Call"]]>0, "AMP", "DELETE")
#   sample_name <- human_cnv[i,][["case"]]
#   each_row <- list(case_name = sample_name,
#                    mut_type = cn_status,
#                    gene_name = gene_list)
#   
#   all_sample_cna <- rbindlist(list(all_sample_cna,each_row))
# }
# all_sample_cna <- setDT(all_sample_cna)
# all_sample_cna <- all_sample_cna[,gene_mutation:=paste(gene_name,mut_type,sep = "_")]
# 
# ### count CNA
# total_cna_gene <- unique(all_sample_cna$gene_mutation)
# 
# all_sample_cna_sum <- list()
# for (each_gene in total_cna_gene){
#   each_sample_number <- nrow(unique(all_sample_cna[gene_mutation == each_gene,.(case_name,gene_mutation)])) 
#   ratio = each_sample_number/human_total_sample_number
#   each_gene_sum <- list(gene_name = each_gene,
#                         numberOfsample = each_sample_number,
#                         ratio = ratio)
#   
#   all_sample_cna_sum <-rbindlist(list(all_sample_cna_sum, each_gene_sum))
#   
# }
# all_sample_cna_sum <- setDT(all_sample_cna_sum)
# all_sample_cna_sum <- all_sample_cna_sum[order(numberOfsample,decreasing = T)]
# all_sample_snv_sum <- all_sample_snv_sum[order(numberOfsample,decreasing = T)]

fwrite(all_sample_snv_sum, file= paste(base_dir,"OM_snv_gene_mutation_sum_03_26.txt",sep =seperator),
       col.names = T, row.names = F, sep = "\t", quote = F, eol="\n")


human_snv_info <- unique(human_snv[,.(Sample,SYMBOL)])
colnames(human_snv_info) <- c("sample_names","gene_name")

human_total_mut <- unique(human_snv_info)
human_total_mut_sample <- unique(human_total_mut$sample_names)

# create a total pathway list (0 and 1 in the list for each sample)
human_total_sum <- list()
for (i in path_col){
  human_total_sum[[i]] =numeric(length(human_total_mut_sample))  
}
human_total_sum[["sample_names"]] =  character(length(human_total_mut_sample))

for (index in 1:length(human_total_mut_sample)){
  each_sample = human_total_mut_sample[index]
  human_total_sum[["sample_names"]][index] = each_sample
  each_sample_gene <- human_total_mut[sample_names==each_sample]$gene_name
  
  for ( col_index in 1:length(path_col)){
    col_name <- path_col[col_index]
    each_path_way_all_gene <- pathway[, col_name, with = F][[col_name]]
    each_path_way_gene <- each_path_way_all_gene[!is.na(each_path_way_all_gene)]
    check_inside <- sum(each_sample_gene %in% each_path_way_gene)
    if (check_inside >0){
      human_total_sum[[col_name]][index] = 1
    }
    else{
      human_total_sum[[col_name]][index] = 0
    }
  }
}

human_total_sum<- setDT(human_total_sum)

## select target tumor type
## human_total_sum is pathway, total_mut is SNV+CNV
human_target_tumor_pathway <- human_total_sum
human_pathway_sum <- NULL
human_total_sample_number <- nrow(human_target_tumor_pathway)
for (each_pathway in colnames(target_pathway)){
  each_pathway_number <- sum(human_target_tumor_pathway[,each_pathway, with =F ][[each_pathway]])
  each_sum <- list(Sample = each_pathway,
                   Mut_type = "path",
                   tumorxsAfflicted = each_pathway_number,
                   tumorsTotal = human_total_sample_number)
  
  human_pathway_sum <- rbindlist(list(human_pathway_sum,each_sum))
}
human_pathway_sum <- setDT(human_pathway_sum)

merge_human_dog <- merge(pathway_sum,human_pathway_sum, by = "Sample")
merge_human_dog <- merge_human_dog[,-5]
colnames(merge_human_dog) <- c("Sample","Mut_type","tumorxsAfflicted","tumorsTotal","humanTumorsAfflicted","humanTumorTotal")


## # ### exam each gene in the target pathway for human
cut_off <- 0
human_target_total_mut <- human_total_mut
human_total_gene_sum <- NULL
for (gene_index in 1:length(target_pathway_gene)){
  each_gene <- target_pathway_gene[gene_index]
  each_gene_tumor_number <- length(unique(human_target_total_mut[gene_name==each_gene]$sample_names))
  if (each_gene_tumor_number >=cut_off){
    each_sum <- list(Sample = each_gene,
                     Mut_type = "inactivating",
                     tumorxsAfflicted = each_gene_tumor_number,
                     tumorsTotal = human_total_sample_number)
    human_total_gene_sum <- rbindlist(list(human_total_gene_sum,each_sum))
  }
  
}
human_gene_list <- human_total_gene_sum$Sample

## # ### exam each gene in the target pathway for dog
cut_off <- 0
dog_total_sample_number <- nrow(dog_target_tumor_pathway)
dog_target_total_mut <- total_mut[Subtype==target_tumor_type,]
dog_total_gene_sum <- NULL
for (gene_index in 1:length(target_pathway_gene)){
  each_gene <- target_pathway_gene[gene_index]
  each_gene_tumor_number <- length(unique(dog_target_total_mut[gene_name==each_gene]$sample_names))
  if (each_gene_tumor_number >=cut_off){
    each_sum <- list(Sample = each_gene,
                     Mut_type = "inactivating",
                     tumorxsAfflicted = each_gene_tumor_number,
                     tumorsTotal = dog_total_sample_number)
    dog_total_gene_sum <- rbindlist(list(dog_total_gene_sum,each_sum))
  }
}
dog_gene_list <- dog_total_gene_sum$Sample
both_dog_human_gene_list <- union(dog_gene_list,human_gene_list)

### examine both human and dog in target_gene list##
cut_off <- 0
human_target_total_mut <- human_total_mut
dog_human_total_gene_sum <- NULL
for (gene_index in 1:length(both_dog_human_gene_list)){
  each_gene <- both_dog_human_gene_list[gene_index]
  human_each_gene_tumor_number <- length(unique(human_target_total_mut[gene_name==each_gene]$sample_names))
  dog_each_gene_tumor_number <- length(unique(dog_target_total_mut[gene_name==each_gene]$sample_names))
  
  each_sum <- list(Sample = each_gene,
                   Mut_type = "inactivating",
                   tumorxsAfflicted = dog_each_gene_tumor_number,
                   tumorsTotal = dog_total_sample_number,
                   humanTumorsAfflicted= human_each_gene_tumor_number,
                   humanTumorTotal = human_total_sample_number)
  dog_human_total_gene_sum <- rbindlist(list(dog_human_total_gene_sum,each_sum))
  
  
}

merge_human_dog <- rbind(merge_human_dog,dog_human_total_gene_sum)
total_row <- nrow(merge_human_dog)
fisher_pval <- c()
for (i in 1:total_row){
  tbl <- matrix(as.numeric(merge_human_dog[i,3:6]), nrow = 2, ncol = 2)
  res <- fisher.test(tbl, alternative = "two.sided")$p.value
  fisher_pval <- c(fisher_pval, res) 
}
merge_human_dog$fisher_pval <- fisher_pval
merge_human_dog <- setDT(merge_human_dog)
merge_human_dog <- merge_human_dog[order(fisher_pval)]

BHcorrect <- p.adjust(merge_human_dog$fisher_pval, method = "BH")
merge_human_dog$BHcorrect<- BHcorrect
merge_human_dog$DogProportion <- merge_human_dog$tumorxsAfflicted/merge_human_dog$tumorsTotal
merge_human_dog$HumanProportion <- merge_human_dog$humanTumorsAfflicted/merge_human_dog$humanTumorTotal
adj_pval <- ifelse(BHcorrect<0.05, "< 0.05", ">=0.05" )
merge_human_dog$Fisher_adj_pval <- adj_pval


fwrite(merge_human_dog, paste(human_base_dir,"Final_human_dog_snv_only_tumor_summary_OM_03_26.txt",sep = seperator),
       quote = F, row.names = F,sep = "\t")

