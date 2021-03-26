source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
## create a target_gene_pathway_list
seperator <- "/"
target_tumor_type="BCL"
human_base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_Cancer_paper/Human/DLBCL"

pathway <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis/PathwayGeneList-03_19.xlsx",
                      sheet = "Pathway_gene_list")
pathway <- setDT(pathway)
target <- c("PI3K","BCR signaling","NFKB","Chromatin remodeler","Cell cycle")
target_pathway <- pathway[,target, with= F]

path_col <- colnames(target_pathway)


target_pathway_gene <- NULL

for (each_pathway in colnames(target_pathway)){
  each_path_gene <- target_pathway[, each_pathway, with = F][[each_pathway]]
  each_path_clean_gene <- each_path_gene[!is.na(each_path_gene)]
  target_pathway_gene <- c(target_pathway_gene,each_path_clean_gene)
}
### analyze dog data

### Dog data ###

base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
seperator <- "/"

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 
#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt")

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))


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

#1. find the most frequeny mut gene for snv and cna in human 
human_total_sample_number <- 42 # the number is from paper
base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_Cancer_paper/Human/TCL"
human_snv <- read.table("clipboard",sep = "\t",header = T,stringsAsFactors = F)
human_snv <- setDT(human_snv)
total_snv_gene <- unique(human_snv$Gene)

## count SNV
all_sample_snv_sum <- list()
for (each_gene in total_snv_gene){
  each_sample_number <- nrow(unique(human_snv[Gene == each_gene,.(Sample,Gene)])) 
  ratio = each_sample_number/human_total_sample_number
  each_gene_sum <- list(gene_name = each_gene,
                        numberOfsample = each_sample_number,
                        ratio = ratio)
  all_sample_snv_sum <-rbindlist(list(all_sample_snv_sum, each_gene_sum))
  
}
all_sample_snv_sum <- setDT(all_sample_snv_sum)

# count snv end 
## prepare CNV sample data
human_cnv <- read.table("clipboard",sep = "\t",header = T,stringsAsFactors = F)
human_cnv <- setDT(human_cnv)
all_sample_cna <- list()
# recreate Data format for cnv
for (i in 1: nrow(human_cnv)){
  gene_list <- unlist(strsplit(human_cnv[i,][['Gene']],","))
  cn_status <- ifelse(human_cnv[i,][["Call"]]>0, "AMP", "DELETE")
  sample_name <- human_cnv[i,][["case"]]
  each_row <- list(case_name = sample_name,
                   mut_type = cn_status,
                   gene_name = gene_list)
  
  all_sample_cna <- rbindlist(list(all_sample_cna,each_row))
}
all_sample_cna <- setDT(all_sample_cna)
all_sample_cna <- all_sample_cna[,gene_mutation:=paste(gene_name,mut_type,sep = "_")]

### count CNA
total_cna_gene <- unique(all_sample_cna$gene_mutation)

all_sample_cna_sum <- list()
for (each_gene in total_cna_gene){
  each_sample_number <- nrow(unique(all_sample_cna[gene_mutation == each_gene,.(case_name,gene_mutation)])) 
  ratio = each_sample_number/human_total_sample_number
  each_gene_sum <- list(gene_name = each_gene,
                        numberOfsample = each_sample_number,
                        ratio = ratio)
  
  all_sample_cna_sum <-rbindlist(list(all_sample_cna_sum, each_gene_sum))
  
}
all_sample_cna_sum <- setDT(all_sample_cna_sum)
all_sample_cna_sum <- all_sample_cna_sum[order(numberOfsample,decreasing = T)]
all_sample_snv_sum <- all_sample_snv_sum[order(numberOfsample,decreasing = T)]

fwrite(all_sample_snv_sum, file= paste(base_dir,"TCL_snv_gene_mutation_sum_03_24.txt",sep =seperator),
       col.names = T, row.names = F, sep = "\t", quote = F, eol="\n")


fwrite(all_sample_cna_sum, file= paste(base_dir,"TCL_cna_gene_mutation_sum_03_24.txt",sep =seperator),
       col.names = T, row.names = F, sep = "\t", quote = F, eol="\n")
################### find the gene mutation frequency end ################### 

################### analyze human data in pathway level ###################   
  
## decide the pathway want to analyze, pathway gene used from oncoprint table S3
## combine snv cnv
human_snv_info <- unique(human_snv[,.(Sample,Gene)])
colnames(human_snv_info) <- c("sample_names","gene_name")
all_sample_cna_info <- unique(all_sample_cna[, .(case_name,gene_name)])
colnames(all_sample_cna_info) <- c("sample_names","gene_name")
human_total_mut <- unique(rbind(human_snv_info,all_sample_cna_info))
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


## dog only
final_dog_only <- rbindlist(list(pathway_sum,dog_total_gene_sum))
final_dog_only <- setDT(final_dog_only)
fwrite(final_dog_only, paste(human_base_dir,"Final_dog_tumor_summary_03_25.txt",sep = seperator),
       quote = F, row.names = F,sep = "\t")

############## Scatter plot ##############
# get the data from human and dog

xlabel <- 'Human CTCL, 42'
ylabel <- 'Canine TCL, 38'

#read.table(paste(pathway_dir,"human_dog_tumor_summary.txt",sep = seperator))

#read.table(paste(pathway_dir,"human_dog_tumor_summary.txt",sep = seperator))
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
pdf(paste(human_base_dir,"humanDogTCL.pdf",sep ="/")
    , height=4.5, width=7.5);

HumanDogPlot<- 
  ggplot(merge_human_dog, 
         aes(x=HumanProportion, y=DogProportion)) +   
  geom_point(aes(size = Fisher_adj_pval 
  )) + 
  ## BH adjusted p-value from fisher test
  scale_color_manual(values=c("black")) +
  scale_size_discrete("Fisher_adj_pval", range=c(1,3)) +
  geom_abline(intercept = 0, color='grey') +
  xlab(xlabel) +
  ylab(ylabel) +
  #   geom_text(data=subset(dat2, Mut_type == "Path"),
  #             aes(label=Gene_pathway ), hjust = 0, nudge_x = 0.02, size = 3) +
  #   geom_text(data=subset(dat2, Human_mut_fraction >= 0.1 & Dog_mut_fraction >= 0.1
  #                         & abs(Human_mut_fraction-Dog_mut_fraction) <= 0.05),  
  #             ## labels all the pathways and genes (mut fraction>0.1, diff<0.05)
  #             aes(label=Gene_pathway), vjust = 0, nudge_x = 0.02, nudge_y = 0.02, check_overlap = TRUE, size = 4) +
  geom_text(data=subset(merge_human_dog, Fisher_adj_pval <= 0.05),
            ## labels all the pathways and genes (mut fraction>0.1, p-value significant)
            aes(label=Sample), vjust = 0, nudge_x = 0.02, nudge_y = 0.02, check_overlap = F, size = 4) +
  scale_y_continuous(breaks = c(0, 0.4, 0.8), limits = c(0, 0.8)) +
  scale_x_continuous(breaks = c(0, 0.4, 0.8), limits = c(0, 0.8)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  theme(legend.position = 'right',
        legend.title = element_text(color = "black", size = 12),
        legend.text = element_text(color = "black", size = 12)) +
  labs(size="Adjusted p-value", colour="Mutation/CNA")
print(HumanDogPlot)
dev.off()

fwrite(merge_human_dog, paste(human_base_dir,"Final_human_dog_tumor_summary_03_24.txt",sep = seperator),
       quote = F, row.names = F,sep = "\t")



