# the script is for table and figure4 in pan-cancer and search Gene from c bioportal not by publication data 
library(data.table)
library(tidyverse)
library(readxl)
library(gridExtra)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
#"/Volumes/Research/GitHub/R_util/my_util.R")
target_tumor_type <- ""
human_total_sample <- 65
## Dog Data ###

base_dir <- 
  #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Oncoprint_analysis"

pathway_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_Cancer_paper/Human/TCL"

seperator <- "/"

whole_wes_clean_breed_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 
#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt")

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))
exclude <- toupper(exclude)
# debug
a = whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID, The_reason_to_exclude, DiseaseAcronym2)]

## pathway data and pathway_gene
pathway <- fread(paste(pathway_dir,"pathway.txt",sep = seperator), na.strings = "")
path_col <- colnames(pathway)

target_pathway_gene <- NULL

for (each_pathway in colnames(pathway)){
  each_path_gene <- pathway[, each_pathway, with = F][[each_pathway]]
  each_path_clean_gene <- each_path_gene[!is.na(each_path_gene)]
  target_pathway_gene <- c(target_pathway_gene,each_path_clean_gene)
}

## add human 5% or >5 mutation gene to examine but only in the target_pathway_gene_list
human_gene_list <- fread(paste(pathway_dir,"Human_gene_list.txt",sep = seperator),header = F, na.strings = "")
human_gene_in_list <- intersect(human_gene_list$V1,target_pathway_gene)

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
total_mut$sample_names <- toupper(total_mut$sample_names)

total_sample <- unique(total_mut$sample_names)
#BCL = unique(whole_wes_clean_breed_table[DiseaseAcronym2=="BCL"]$Case_ID)
total_mut <- total_mut[!sample_names %in% exclude,]
#tot_mut_bcl <- unique(total_mut[Subtype=="BCL"]$sample_names)
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
dog_target_tumor_pathway <- total_sum[Subtype==target_tumor_type & !sample_names %in% exclude,]
pathway_sum <- NULL
total_tumor_number <- nrow(dog_target_tumor_pathway)
for (each_pathway in colnames(pathway)){
  each_pathway_number <- sum(dog_target_tumor_pathway[,each_pathway, with =F ][[each_pathway]])
  each_sum <- list(Sample = each_pathway,
                   Mut_type = "path",
                   tumorxsAfflicted = each_pathway_number,
                   tumorsTotal = total_tumor_number)
                    
  pathway_sum <- rbindlist(list(pathway_sum,each_sum))
}
pathway_sum <- setDT(pathway_sum)
# 
# ### exam each gene in the target pathway
cut_off <- 5
target_total_mut <- total_mut[Subtype==target_tumor_type,]
total_gene_sum <- NULL
for (gene_index in 1:length(target_pathway_gene)){
  each_gene <- target_pathway_gene[gene_index]
  each_gene_tumor_number <- length(unique(target_total_mut[gene_name==each_gene]$sample_names))
  if (each_gene_tumor_number >=cut_off){
  each_sum <- list(Sample = each_gene,
                   Mut_type = "inactivating",
                   tumorxsAfflicted = each_gene_tumor_number,
                   tumorsTotal = total_tumor_number)
  total_gene_sum <- rbindlist(list(total_gene_sum,each_sum))
  }
  else{
    if (each_gene %in% human_gene_in_list){
      each_sum <- list(Sample = each_gene,
                       Mut_type = "inactivating",
                       tumorxsAfflicted = each_gene_tumor_number,
                       tumorsTotal = total_tumor_number)
      total_gene_sum <- rbindlist(list(total_gene_sum,each_sum))
    }
  }
}

final_sum <- rbindlist(list(pathway_sum,total_gene_sum))
final_sum <- setDT(final_sum)

fwrite(final_sum, paste(pathway_dir,"dog_tumor_pathway_summary.txt",sep = seperator),
       quote = F, row.names = F,sep = "\t")


###### Human data  if cbioprotal has, then I use the number in there ## ####
# human_base <- "G:/MAC_Research_Data/Pan_cancer/Pan_Cancer_paper/Human/DLBCL/dlbc_tcga_pan_can_atlas_2018"
# seperator <- "/"
# human_data <- fread(paste(pathway_dir,"gene_alteration_frequency.tsv",sep = seperator))
# human_total_cnv <- fread(paste(human_base,"data_CNA.txt",sep = seperator))
# human_total_sample <- colnames(human_total_cnv)[-c(1,2)]



# 
# target_human <- human_data[,.(`Gene Symbol`,`Num Samples Altered`)]
# colnames(target_human) <- c("Sample","humanTumorsAfflicted")
# target_human <- data.table(target_human,key = "Sample")
# #merge(final_sum,target_human, by.x = "Sample", by.y="Gene Symbol")
# final_output <- target_human[final_sum]
# final_output <- final_output[,humanTumorTotal:=human_total_sample]
# setcolorder(final_output,c("Sample","Mut_type","tumorxsAfflicted","tumorsTotal","humanTumorsAfflicted","humanTumorTotal"))
# fwrite(final_output,file = paste(pathway_dir,"human_dog_tumor_summary.txt",sep = seperator),
#        col.names = T, row.names = F, quote = F,sep = "\t", na = "NA")



############## Scatter plot ##############
# get the data from human and dog

xlabel <- 'Human DLBCL,48'
ylabel <- 'Canine BCL, 55'
HumanDogPathway <- read.table("clipboard",sep = "\t",header = T,stringsAsFactors = F)
  #read.table(paste(pathway_dir,"human_dog_tumor_summary.txt",sep = seperator))
total_row <- nrow(HumanDogPathway)
fisher_pval <- c()
for (i in 1:total_row){
  tbl <- matrix(as.numeric(HumanDogPathway[i,3:6]), nrow = 2, ncol = 2)
  res <- fisher.test(tbl, alternative = "two.sided")$p.value
  fisher_pval <- c(fisher_pval, res) 
}
HumanDogPathway$fisher_pval <- fisher_pval
HumanDogPathway <- setDT(HumanDogPathway)
HumanDogPathway <- HumanDogPathway[order(fisher_pval)]

BHcorrect <- p.adjust(HumanDogPathway$fisher_pval, method = "BH")
HumanDogPathway$BHcorrect<- BHcorrect
HumanDogPathway$DogProportion <- HumanDogPathway$tumorxsAfflicted/HumanDogPathway$tumorsTotal
HumanDogPathway$HumanProportion <- HumanDogPathway$humanTumorsAfflicted/HumanDogPathway$humanTumorTotal
adj_pval <- ifelse(BHcorrect<0.05, "< 0.05", ">=0.05" )
HumanDogPathway$Fisher_adj_pval <- adj_pval
pdf(paste(pathway_dir,"humanDogBCL.pdf",sep ="/")
    , height=4.5, width=7.5);



HumanDogPlot<- 
  ggplot(HumanDogPathway, 
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
    geom_text(data=subset(HumanDogPathway, Fisher_adj_pval <= 0.01),
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

fwrite(HumanDogPathway, paste(pathway_dir,"Final_human_dog_tumor_summary.txt",sep = seperator),
          quote = F, row.names = F,sep = "\t")

### all together
total_data<- read.table("clipboard",sep = "\t",header = T,stringsAsFactors = F)
total_data <- setDT(total_data)
pdf(paste(pathway_dir,"humanDogall.pdf",sep ="/")
    , height=9, width=16);

x_total_label <- c('Human DLBCL,48',"Human TCL, 42",'Human, 65')
y_total_label <- c('Canine BCL, 55',"Canine TCL, 38",'Canine 71')
tumor_type <- c("BCL","TCL","OM")


plot <- list()
for (i in 1:length(tumor_type)){
  xlabel <- x_total_label[i]
  ylabel <- y_total_label[i]
  each_tumor_type <- tumor_type[i]
  HumanDogPathway <-   total_data[Tumor==each_tumor_type,]
  
HumanDogPlot<- 
  ggplot(HumanDogPathway, 
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
  geom_text(data=subset(HumanDogPathway, Fisher_adj_pval <= 0.01),
            ## labels all the pathways and genes (mut fraction>0.1, p-value significant)
            aes(label=Pathway), vjust = 0, nudge_x = 0.02, nudge_y = 0.02, check_overlap = F, size = 4) +
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

plot[[i]] <- HumanDogPlot
#print(HumanDogPlot)
}
length(plot)
do.call("grid.arrange", c(plot, ncol=2))
dev.off()

