library(tidyverse)
library(data.table)
library(gridExtra)
library(ggrepel)
library(readxl)
library(ggpubr)
library(grid)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
  #"/Volumes/Research/GitHub/R_util/my_util.R")

source("C:/Users/abc73/Documents/GitHub/VAF/six_base_function_util.R")
#source("/Volumes/Research/GitHub/VAF/six_base_function_util.R")
seperator <- "/"
whole_wes_clean_breed_table <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt")
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt") 

fontsize=35
dot_size <- 1.4;
abs_text_size <- 50;
regular.text <- element_text(colour="black",size=abs_text_size);
xangle <- 45;
xaxis_just <- ifelse(xangle > 0, 1, 0.5);

### add extra parameter (use 1. a+b-a&b 2. min of a or b, 3. average a or b)
create_overlap_summary <- function(our_data,publis_data,intercet_sample, method ){
  
  total_uniq_num_to_them <- NULL
  total_uniq_num_to_us <- NULL
  total_share_number <- NULL
  total_uniq_ratio_to_them <- NULL
  total_uniq_ratio_to_us <- NULL
  total_share_ratio <- NULL
  total_sample <- NULL
  total_denomitor <- NULL
  total_their_mut_number <- NULL
  total_UGA_mut_number <- NULL
  for (samp in intercet_sample){
    
    our_each_mut <- our_data[Case == samp, .(chrom_loc)]
    publish_each_mut <- publis_data[Case == samp,.(chrom_loc)]
    
    
    our_each <- nrow(our_each_mut)
    sanger_each <- nrow(publish_each_mut)
    intercet_data <- nrow(intersect(our_each_mut,publish_each_mut))
    
    if (method == 'min'){
      denominator <- min(c(our_each,sanger_each))
    }
    
    else if (method == 'average'){
      denominator <- mean(c(our_each,sanger_each))
    }
    else if (method == 'union'){
      denominator <- our_each+sanger_each-intercet_data
    }
    
    # count
    number_overlap <- nrow(intersect(our_each_mut,publish_each_mut))
    uniq_number_to_us <-nrow(setdiff(our_each_mut,publish_each_mut)) 
    uniq_number_to_them <- nrow(setdiff(publish_each_mut,our_each_mut)) 
    their_mut_number <- nrow(publish_each_mut)
    UGA_mut_number <- nrow(our_each_mut)
    # ratio
    uniq_ratio_to_them <- uniq_number_to_them/(denominator)
    uniq_ratio_to_us <- uniq_number_to_us/(denominator)
    overlap_ratio <- number_overlap/(denominator)
    
    
    # summary
    total_share_ratio <- c(total_share_ratio,overlap_ratio)
    total_uniq_ratio_to_us <- c(total_uniq_ratio_to_us,uniq_ratio_to_us)
    total_uniq_ratio_to_them <- c(total_uniq_ratio_to_them,uniq_ratio_to_them)
    total_uniq_num_to_them <- c(total_uniq_num_to_them,uniq_number_to_them)
    total_uniq_num_to_us <- c(total_uniq_num_to_us,uniq_number_to_us)
    total_share_number <- c(total_share_number,number_overlap)
    total_sample <- c(total_sample,samp)
    total_denomitor <- c(total_denomitor,denominator)
    total_their_mut_number <- c(total_their_mut_number,their_mut_number)
    total_UGA_mut_number <- c(total_UGA_mut_number,UGA_mut_number)
    
  }
  
  data <- data.frame(share_ratio = as.numeric(total_share_ratio), 
                     sample = total_sample,
                     uniq_ratio_to_uga = as.numeric(total_uniq_ratio_to_us),
                     uniq_ratio_to_publication = as.numeric(total_uniq_ratio_to_them),
                     Unique_to_original_publication = as.numeric(total_uniq_num_to_them),
                     Unique_to_our_study = as.numeric(total_uniq_num_to_us),
                     Shared = as.numeric(total_share_number),
                     total_denomitor= as.numeric(total_denomitor),
                     total_their_mut_number=as.numeric(total_their_mut_number),
                     total_UGA_mut_number=as.numeric(total_UGA_mut_number))
  data <- setDT(data)
  
  return (data)
}


create_overlap_summary_for_each <- function(our_data,publis_data,intercet_sample){
  final_sum <- list()
  unique_to_them_sum <- data.table()
  unique_to_us_sum <- data.table()
  share_sum <- data.table()
  for (samp in intercet_sample){
    our_each_mut <- unique(our_data[Case == samp, .(chrom_loc)][['chrom_loc']])
    publish_each_mut <- unique(publis_data[Case == samp,.(chrom_loc)][['chrom_loc']])
    intercet_data <- intersect(our_each_mut,publish_each_mut)
    each_unique_to_us <- setdiff(our_each_mut,publish_each_mut)
    each_unique_to_publish <- setdiff(publish_each_mut,our_each_mut)
    
    if (length(each_unique_to_publish)==0){
      each_unique_to_publish <- "NA"
    }
    if (length(each_unique_to_us)==0){
      each_unique_to_us <- "NA"
    }
    if (length(intercet_data)==0){
      intercet_data <- "NA"
    }
    each_unique_to_us <- list(Case = samp,
                              unique_to_us =each_unique_to_us)
    
    each_unique_to_publish <- list(Case = samp,
                                   unique_to_publish =each_unique_to_publish)
    each_share <- list(Case = samp,
                       share = intercet_data)
    
    unique_to_them_sum <- rbindlist(list(unique_to_them_sum,each_unique_to_publish))
    
    unique_to_us_sum <- rbindlist(list(unique_to_us_sum,each_unique_to_us))
    
    share_sum <- rbindlist(list(share_sum,each_share))
    
    # summary
    
  }
  unique_to_them_sum <- setDT(unique_to_them_sum)
  unique_to_us_sum <- setDT(unique_to_us_sum)
  share_sum <- setDT(share_sum)
  final_sum[['unique_to_them_sum']] <- unique_to_them_sum
  final_sum[['unique_to_us_sum']] <- unique_to_us_sum
  final_sum[['share_sum']] <- share_sum
  return (final_sum)
}



##### check overlap samples ######
base <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Compare_publication/OM_mutation_compare_with_Sanger"
  "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Compare_publication\\OM_mutation_compare_with_Sanger"

## clean sanger_data
original_signature <- read_excel(paste(base,"Sanger_mutation.xlsx",sep =seperator),
                                 sheet ="Canine", skip = 29)
original_signature <- setDT(original_signature)
original_signature$Chromosome <- paste("chr",original_signature$`#Chr`,sep="")
original_signature$chrom_loc <- paste(original_signature$Chromosome,original_signature$Position,sep="_")
original_signature$merge_key <- paste(original_signature$Chromosome,original_signature$Position,original_signature$Ref,original_signature$Alt,sep="_")
target_original <- original_signature[,.(Sample,Chromosome,Position,Ref,Alt,Consequence,chrom_loc,merge_key)]
#unique(target_original$Consequence)
#target_conseq <- unique(target_original$Consequence[!grepl("frame",target_original$Consequence)])

seperator <- "/"
sanger_signature <- fread(paste(base,"DbSNP_sanger_CDS_mut_file_after_DBSNP_03_23.txt",sep =seperator))
sanger_filter_SNP <- sanger_signature
sanger_signature[Sample=="DD0013a",]

sanger_filter_SNP$chrom_loc <- paste(sanger_filter_SNP$Chromosome,sanger_filter_SNP$Position,sep = "_")
sanger_filter_SNP$merge_key <- paste(sanger_filter_SNP$Chromosome,sanger_filter_SNP$Position,sanger_filter_SNP$Ref,sanger_filter_SNP$Alt,sep="_")

samples <- sort(unique(sanger_filter_SNP$Sample))

sanger_filter_SNP[Sample=="DD0013a",]
merge_table <-  merge(x = sanger_filter_SNP, y = target_original, by.x="merge_key", by.y="merge_key", all.x = T)
target_col <- c("Sample.x","Chromosome.x","Position.x","Ref.x","Alt.x","chrom_loc.x","Consequence")
merge_table <- merge_table[,target_col,with = F]
a = merge_table[duplicated(merge_table),]
merge_table <-  distinct(merge_table)
colnames(merge_table) <- c("Sample","Chromosome","Position","Ref","Alt","chrom_loc","Consequence")
## only select SNV data ###
#target_conseq <- c("missense_variant","synonymous_variant","stop_gained","start_lost")

#unique(target_original$Consequence)

final_merge <- merge_table
clean_sanger <- setDT(final_merge)
colnames(clean_sanger)[1] <- "Case"
# fwrite(clean_sanger, file = paste(base,"DbSNP_sanger_snv_indel_CDS_mut_file_03_30.txt",sep = seperator),
#        col.names = T, row.names = F, sep = "\t",quote = F,eol = "\n")
#### Analyzed the before 5steps

# ## indel file
indel_file <- fread(paste(base,"total_CDS_indel_info_withGene_04_08.txt",sep =seperator))
indel_col_names <- c("chrom","pos","ref","alt","gene_name","ensembl_id","status","sample_names")
colnames(indel_file) <- indel_col_names
indel_file$tumor_type <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table)
final_indel <- indel_file[tumor_type=="OM",.(sample_names,chrom,pos,ref,alt,status)]
final_indel$chrom_loc <- paste(final_indel$chrom,final_indel$pos,sep ="_")
final_indel$Case <-  sapply(final_indel$sample_names,convert_sample)



our_OM_before <- fread(paste(base,"DbSNP_CDS_Total_Before_5step_Sanger_Mutect1_03_29.txt",sep = seperator));
# colnames(our_OM_before) <- c("chrom","pos","vaf","ref","alt","tRef","tAlt","sample_names","gene_name","ensembl_id",
#                "status","tumor_type","symbol")
our_OM_before <- our_OM_before[symbol=="OM SC",]
our_OM_before <- our_OM_before[tumor_type=="OM",.(sample_names,chrom,pos,ref,alt,status)]
our_OM_before$Case <- sapply(our_OM_before$sample_names,convert_sample)
our_OM_before$chrom_loc <- paste(our_OM_before$chrom,our_OM_before$pos,sep ="_")
our_OM_before <- rbind(our_OM_before,final_indel)
our_sample <- sort(unique(our_OM_before$Case))
their_sample <- sort(unique(clean_sanger$Case))
OM_before_intercet_sample <- intersect(their_sample,our_sample)
fwrite(our_OM_before, file = paste(base, "UGA_Mutect1_before_steps_withStrelka_03_31.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t",eol = "\n")

### after 5 steps
our_OM_after <- fread(paste(base,"DbSNP_CDS_Total_After_5step_Sanger_Mutect1_03_29.txt",sep = seperator));
# colnames(our_OM_after) <- c("chrom","pos","vaf","ref","alt","tRef","tAlt","sample_names","gene_name","ensembl_id",
#                              "status","tumor_type","symbol")

our_OM_after[Case=='DD0041a',]$chrom_loc

our_OM_after <- our_OM_after[symbol=="OM SC",.(sample_names,chrom,pos,ref,alt,status)]
our_OM_after$Case <- sapply(our_OM_after$sample_names, convert_sample)
our_OM_after$chrom_loc <- paste(our_OM_after$chrom,our_OM_after$pos,sep= "_")
#our_OM_after <- clean_table(our_OM_after)
our_OM_after$chrom_loc <- paste(our_OM_after$chrom,our_OM_after$pos,sep ="_")
our_OM_after <- rbind(our_OM_after,final_indel)
our_sample <- sort(unique(our_OM_after$Case))
their_sample <- sort(unique(clean_sanger$Case))
OM_after_intercet_sample <- intersect(their_sample,our_sample)
fwrite(our_OM_after, file = paste(base, "UGA_Mutect1_after5_steps_withStrelka_03_31.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t",eol = "\n")

## Burair
our_OM_Burair <- fread(paste(base,"Final_Total_withGene_Burair_Filtering3_VAF_Mutect_orientBiasModified_04_02txt.gz",sep = seperator));
our_OM <- our_OM_Burair[symbol=="OM SC",.(sample_names,chrom,pos,ref,alt,status)]

our_OM$Case <- sapply(our_OM$sample_names, convert_sample)
our_OM$chrom_loc <- paste(our_OM$chrom,our_OM$pos,sep= "_")
our_OM <- rbind(our_OM,final_indel)

fwrite(our_OM, file = paste(base, "Final_burair_pipeline_withStrelka_04_13.txt",sep = seperator),
       col.names = T, row.names = F, quote = F, sep = "\t",eol = "\n")

our_sample <- unique(our_OM$Case)
their_sample <- sort(unique(clean_sanger$Case))
Burair_intercet_sample <- intersect(their_sample,our_sample)
Burair_diff <- setdiff(their_sample,our_sample)

total_three_intercet <- intersect(intersect(OM_before_intercet_sample,OM_after_intercet_sample),
                                  Burair_intercet_sample)

##### check overlap samples end ######

# our_OM_before <- fread(paste(base,"total_final_without_Gene_Burair_Filtering3_VAF_Mutect_Before_0201.txt.gz",sep = seperator));
# our_OM_before <- our_OM_before[symbol=="OM SC",]

png(file = paste(base,"04_26","include_indel_5steps_only_Mutation_number_compare_with_OM_publication_03_29.png",sep =seperator),
    width = 5000, height =2700, units = "px", res = 500)

data_5setps <- create_overlap_summary(our_OM_after,clean_sanger,total_three_intercet,method = 'union')
after_5steps_data_sum <- create_overlap_summary_for_each(our_OM_after,clean_sanger,total_three_intercet)

# fwrite(after_5steps_data_sum$unique_to_them_sum, file = paste(base,"unique_to_them.txt",sep = "/"),
#        col.names = T, row.names = F,sep = "\t",eol = "\n")
# 
# fwrite(after_5steps_data_sum$unique_to_us_sum, file = paste(base,"unique_to_us.txt",sep = "/"),
#        col.names = T, row.names = F,sep = "\t",eol = "\n")
# fwrite(after_5steps_data_sum$share_sum, file = paste(base,"share_both.txt",sep = "/"),
#        col.names = T, row.names = F,sep = "\t",eol = "\n")

count_data <- melt(data_5setps, id.vars = c("sample"),
                   measure.vars= c("Unique_to_our_study","Unique_to_original_publication","Shared"),
                   variable.name = "fill")
count_data <- count_data[order(sample)]

x <- count_data$sample
y <- count_data$value
classify <- c("Unique_to_our_study","Unique_to_original_publication","Shared");
fill <- count_data$fill
fill <- factor(fill, levels=classify);
samples <- unique(x);
sample_order <- order(sapply(samples, function(s) {sum(y[which(x == s)])}), decreasing=T);
after_5steps_order <- sample_order
x <- factor(x, levels=samples[sample_order]);
plot_data <- data.frame(x=x, y=y, fill=fill);

fill_colors <- c("cyan","black","red");

p <- my_bar_function(plot_data,fill_colors = fill_colors,
                     title="Somatic base substitutions\n after 5-step filtering",fontsize=35)
p <- p+scale_y_continuous(breaks=c(0,50,100,150))
p <- p+theme(legend.position="none",
             axis.text=regular.text, 
             axis.title=regular.text)
print(p)
fwrite(data_5setps,file=paste(base,"04_26","include_indel_5steps_only_compare_publication.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
dev.off()
# 
# png(file = paste(base,"04_26","5steps_only_Mutation_ratio_compare_with_OM_publication_03_29.png",sep =seperator),
#     width = 5000, height =2700, units = "px", res = 500)
# 
# ratio_data <- melt(data_5setps, id.vars = c("sample"),
#                    measure.vars= c("uniq_ratio_to_uga","uniq_ratio_to_publication","share_ratio"),
#                    variable.name = "fill")
# 
# ratio_data <- ratio_data[order(sample)]
# 
# x <- ratio_data$sample
# y <- ratio_data$value
# classify <- c("uniq_ratio_to_uga","uniq_ratio_to_publication","share_ratio");
# fill <- ratio_data$fill
# fill <- factor(fill, levels=classify);
# samples <- unique(x);
# sample_order <-sample_order;
# x <- factor(x, levels=samples[sample_order]);
# plot_data <- data.frame(x=x, y=y, fill=fill);
# #plot_data <- plot_data[-which(plot_data$x =="CMT-033"),] # remove the outlier
# fill_colors <- c("cyan","black","red");
# 
# p1 <- my_bar_function(plot_data,fill_colors = fill_colors,
#                       title="UGA OM mutation ratio overlapped with Sanger Publication",fontsize=20)
# print(p1)
# dev.off()

### 5steps only end

## Before 5 steps

png(file = paste(base,"04_26","include_indel_Before_5_steps_Mutation_number_compare_with_OM_publication_03_29.png",sep =seperator),
    width = 5000, height =2700, units = "px", res = 500)

data_before_5setps <- create_overlap_summary(our_OM_before,clean_sanger,total_three_intercet,method = 'union')

count_data <- melt(data_before_5setps, id.vars = c("sample"),
                   measure.vars= c("Unique_to_our_study","Unique_to_original_publication","Shared"),
                   variable.name = "fill")
count_data <- count_data[order(sample)]

x <- count_data$sample
y <- count_data$value
classify <- c("Unique_to_our_study","Unique_to_original_publication","Shared");
fill <- count_data$fill
fill <- factor(fill, levels=classify);
samples <- unique(x);
sample_order <- after_5steps_order
x <- factor(x, levels=samples[sample_order]);
plot_data <- data.frame(x=x, y=y, fill=fill);

fill_colors <- c("cyan","black","red");

p <- my_bar_function(plot_data,fill_colors = fill_colors,
                     title="Somatic base substitutions\nvia MuTect alone",fontsize=35)

p <- p+theme(legend.position="none",
             axis.text=regular.text, 
             axis.title=regular.text)
             
print(p)

fwrite(data_before_5setps,file=paste(base,"04_26","include_indel_Before_5steps_Mutect1_compare_publication.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")

dev.off()
# 
# png(file = paste(base,"04_26","Before_5_steps_Mutation_ratio_compare_with_OM_publication_03_29.png",sep =seperator),
#     width = 5000, height =2700, units = "px", res = 500)
# 
# ratio_data <- melt(data_before_5setps, id.vars = c("sample"),
#                    measure.vars= c("uniq_ratio_to_uga","uniq_ratio_to_publication","share_ratio"),
#                    variable.name = "fill")
# 
# ratio_data <- ratio_data[order(sample)]
# 
# x <- ratio_data$sample
# y <- ratio_data$value
# classify <- c("uniq_ratio_to_uga","uniq_ratio_to_publication","share_ratio");
# fill <- ratio_data$fill
# fill <- factor(fill, levels=classify);
# samples <- unique(x);
# sample_order <-sample_order;
# x <- factor(x, levels=samples[sample_order]);
# plot_data <- data.frame(x=x, y=y, fill=fill);
# #plot_data <- plot_data[-which(plot_data$x =="CMT-033"),] # remove the outlier
# fill_colors <- c("cyan","black","red");
# 
# p1 <- my_bar_function(plot_data,fill_colors = fill_colors,
#                       title="UGA OM mutation ratio overlapped with Sanger Publication",fontsize=20)
# print(p1)


dev.off()

####### Burair filtering #########


png(file = paste(base,"04_26","Burair_filtering_Mutect1_Mutation_number_compare_with_OM_publication_03_29.png",sep =seperator),
    width = 5000, height =2700, units = "px", res = 500)

Burair_filtering_data <- create_overlap_summary(our_OM,clean_sanger,total_three_intercet, method = 'union')

count_data <- melt(Burair_filtering_data, id.vars = c("sample"),
                   measure.vars= c("Unique_to_our_study","Unique_to_original_publication","Shared"),
                   variable.name = "fill")
count_data <- count_data[order(sample)]

x <- count_data$sample
y <- count_data$value
classify <- c("Unique_to_our_study","Unique_to_original_publication","Shared");
fill <- count_data$fill
fill <- factor(fill, levels=classify);
samples <- unique(x);
sample_order <- after_5steps_order
x <- factor(x, levels=samples[sample_order]);
plot_data <- data.frame(x=x, y=y, fill=fill);

fill_colors <- c("cyan","black","red");

p <- my_bar_function(plot_data,fill_colors = fill_colors,
                     title="Somatic base substitutions\n after 5-step & orientation bias filtering",
                     fontsize=35)
p <- p+scale_y_continuous(breaks=c(0,50,100,150))
p <- p+theme(legend.position="none",
             axis.text=regular.text, 
             axis.title=regular.text)
print(p)

dev.off()
fwrite(Burair_filtering_data,file=paste(base,"04_26","Burair_filtering_Mutect1_compare_publication.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
# 
# png(file = paste(base,"04_26","Burair_filtering_Mutect1_Mutation_ratio_compare_with_OM_publication_03_29.png",sep =seperator),
#     width = 5000, height =2700, units = "px", res = 500)
# 
# ratio_data <- melt(Burair_filtering_data, id.vars = c("sample"),
#                    measure.vars= c("uniq_ratio_to_uga","uniq_ratio_to_publication","share_ratio"),
#                    variable.name = "fill")
# 
# ratio_data <- ratio_data[order(sample)]
# 
# x <- ratio_data$sample
# y <- ratio_data$value
# classify <- c("uniq_ratio_to_uga","uniq_ratio_to_publication","share_ratio");
# fill <- ratio_data$fill
# fill <- factor(fill, levels=classify);
# samples <- unique(x);
# sample_order <-sample_order;
# x <- factor(x, levels=samples[sample_order]);
# plot_data <- data.frame(x=x, y=y, fill=fill);
# #plot_data <- plot_data[-which(plot_data$x =="CMT-033"),] # remove the outlier
# fill_colors <- c("cyan","black","red");
# 
# p1 <- my_bar_function(plot_data,fill_colors = fill_colors,
#                       title="UGA OM mutation ratio overlapped with Sanger Publication",fontsize=20)
# print(p1)




### compare end ###

# 
# ############# Plot the sanger six bases #############
# 
# 
# fontsize <- 20;
# signature_colors <- c("cyan","black","red","gray","green","pink");
# signature_levels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
# col_name <- c("number","ref","alt")
# OM_base <- "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Six_base_sub_FFPE"
# MT_base <- "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Six_base_sub_FFPE"
# OM_sample <-  sort(list.files (path =paste(OM_base,"OM","Mutect1",sep= "\\")))
# MC_sample <- sort(list.files (path = paste(MT_base,"MT","Mutect1",sep= "\\")))
# fill_colors <- signature_colors 
# Cancer_type <- "OM"
# 
# ## to see if we want to excldued
# # MC_sample <- MC_sample[-match('CMT-33', MC_sample)]
# xangle <- 45
# 
# if (Cancer_type =="OM"){
#   total_sample <-OM_sample
#   base <- OM_base
# }else{
#   total_sample <- MC_sample
#   base <- MT_base
# }
# mut_before_count_sum <- NULL
# mut_after_count_sum <- NULL
# mut2_before_count_sum <- NULL
# mut2_after_count_sum <- NULL
# 
# sanger_signature <- read_excel(paste(base,"Sanger_mutation.xlsx",sep="\\"),
#                                sheet ="Canine", skip = 29)
# sanger_signature <- setDT(sanger_signature)
# sample <- sort(unique(sanger_signature$Sample))
# clean_sanger <- sanger_signature[, .(Ref,Alt,Sample)]
# colnames(clean_sanger) <- c("ref","alt","sample_name")
# 
# clean_sanger <- clean_table(clean_sanger)
# 
# clean_sanger$conver_mut_type <- sapply(clean_sanger$mut_type, convert_mutation_type)
# 
# # pdf(paste(OM_base,"six_base_Compare_OM.pdf",sep="\\")
# #     , height=8.94, width=12.84);
# 
# 
# for (sample in total_sample){
#   if(grepl("-1",sample)){
#     split_words <- str_split(sample,'-')[[1]][1]
#     sanger_sample <- paste(split_words,'c',sep = "")
#   }else if(grepl("-2",sample)){
#     split_words <- str_split(sample,'-')[[1]][1]
#     sanger_sample <- paste(split_words,'d',sep = "")
#   }
#   else{
#     sanger_sample <- paste(sample,'a',sep = "")
#   }
#   
#   ## calculate Sanger mut and 6 bases
#   # data <- as.data.frame(table(clean_sanger[sample_name== sanger_sample,.(conver_mut_type)]))
#   # colnames(data) <- c("conver_mut_type","number")
#   # data <- setDT(data)
#   # data_rate <- count_mutation_rates(data,signature_levels)
#   # data_number <- count_mutation_number(data,signature_levels)
#   # rate_data <- prepare_bar_plot_from_table(data_rate,sanger_sample)
#   # number_data <- prepare_bar_plot_from_table(data_number,sanger_sample)
#   # 
#   # ps <- my_barplot(rate_data, fill_colors=signature_colors, abs_text_size=12, xangle=30)
#   # ps <- ps+labs(title="Sanger_mut_Ratio", y="", x="", fill="Mut.");
#   # 
#   # psn <- my_barplot(number_data, fill_colors=signature_colors, abs_text_size=12, xangle=30)
#   # psn <- psn+labs(title="Sanger_mut_Count", y="", x="", fill="Mut.");
#   
#   ## calculate our own data
#   M1before_file_name <- paste(sample,"Mutect1_before5steps.txt",sep = '_')
#   M1after_file_name <- paste(sample,"Mutect1_after5steps.txt",sep = '_')
#   M2before_file_name <- paste(sample,"Mutect2_before5steps.txt",sep = '_')
#   M2after_file_name <- paste(sample,"Mutect2_after5steps.txt",sep = '_')
#   # #
#   Mutect1_before <-fread(paste(base,Cancer_type,"Mutect1",sample,"before",M1before_file_name ,sep="\\"))
#   Mutect1_after <- fread(paste(base,Cancer_type,"Mutect1",sample, "after",M1after_file_name ,sep="\\"))
#   Mutect2_before <-fread(paste(base,Cancer_type,"Mutect2",sample,"before",M2before_file_name ,sep="\\"))
#   Mutect2_after <- fread(paste(base,Cancer_type,"Mutect2",sample, "after",M2after_file_name ,sep="\\"))
#   # #
#   mut_before_rate <- find_six_base(Mutect1_before,sample,rate =T);
#   mut_after_rate <- find_six_base(Mutect1_after,sample,rate =T);
#   mut2_before_rate <- find_six_base(Mutect2_before,sample,rate =T);
#   mut2_after_rate <- find_six_base(Mutect2_after,sample,rate =T);
#   # # 
#   mut_before_count <- find_six_base(Mutect1_before,rate =F,sample_name=sample);
#   mut_after_count <- find_six_base(Mutect1_after,rate =F,sample_name=sample);
#   mut2_before_count <- find_six_base(Mutect2_before,rate =F,sample_name=sample);
#   mut2_after_count <- find_six_base(Mutect2_after,rate =F,sample_name=sample);
#   # # 
#   mut_before_count_sum <- rbind(mut_before_count_sum,mut_before_count)
#   mut_after_count_sum <- rbind(mut_after_count_sum,mut_after_count)
#   mut2_before_count_sum <- rbind(mut2_before_count_sum,mut2_before_count)
#   mut2_after_count_sum <- rbind(mut2_after_count_sum,mut2_after_count)
#   # 
#   # 
#   # ## Plot the signature
#     # p1 <- my_barplot(mut_before_rate, fill_colors, abs_text_size=12, xangle=30)
#     # p1 <- p1 + labs(title="mut_before_Ratio", y="", x="", fill="Mut.");
#     # p2 <- my_barplot(mut_after_rate,fill_colors=signature_colors, abs_text_size=12, xangle=30)
#     # p2 <- p2 + labs(title="mut_after_Ratio", y="", x="", fill="Mut.");
#     # p3 <- my_barplot(mut2_before_rate,fill_colors, abs_text_size=12, xangle=30)
#     # p3 <- p3 + labs(title="mut2_before_Ratio", y="", x="", fill="Mut.");
#     # p4 <- my_barplot(mut2_after_rate,fill_colors=signature_colors, abs_text_size=12, xangle=30)
#     # p4 <- p4 + labs(title="mut2_after_Ratio", y="", x="", fill="Mut.");
#   # # #   
#   # # #   
#     # p1_1 <- my_barplot(mut_before_count,fill_colors, abs_text_size=12, xangle=30)
#     # p1_1 <- p1_1 + labs(title="mut_before_Counts", y="", x="", fill="Mut.");
#     # p2_1 <- my_barplot(mut_after_count,fill_colors=signature_colors, abs_text_size=12, xangle=30)
#     # p2_1 <- p2_1 + labs(title="mut_after_Counts", y="", x="", fill="Mut.");
#     # p3_1 <- my_barplot(mut2_before_count,fill_colors, abs_text_size=12, xangle=30)
#     # p3_1 <- p3_1 + labs(title="mut2_before_Counts", y="", x="", fill="Mut.");
#     # p4_1 <- my_barplot(mut2_after_count,fill_colors=signature_colors, abs_text_size=12, xangle=30)
#     # p4_1 <- p4_1 + labs(title="mut2_after_Counts", y="", x="", fill="Mut.");
# 
#   # # #   lay <- rbind(c(1,2),
#   # # #                c(3,4))
#   # # #   #c(6,7,8,9,9))
#     # grid.arrange(p1,p2,p1_1,p2_1,
#     #              p3,p4,p3_1,p4_1,
#     #              nrow = 2, top = textGrob(sample,gp=gpar(fontsize=24,font=3)))
# 
# }
# 
# dev.off()
# 
# ## check what samples are not normal
# ex <- setDT(mut_before_count_sum)
# answer <- ex[, .(total=sum(y)),keyby = .(x)][order(-total)]
# 
# ### check the sample order
# # mut_before_count_sum 
# # mut_after_count_sum 
# # mut2_before_count_sum 
# # mut2_after_count_sum 
# 
# my_new_bar_plot <- function(data,fill_colors,title){
#   p <- ggplot(data, aes(x=x, y=y, fill=fill)) + geom_bar(stat="identity",position='stack', width=0.6)+
#     ggtitle(title)+ 
#     scale_fill_manual(values=fill_colors)+
#     theme(
#       plot.title = element_text(size = 20, face = "bold"),
#       axis.title.x = element_blank(),
#       #element_text(face="plain",colour="black",size=fontsize),
#       axis.title.y = element_blank(),
#       #element_text(face="plain",colour="black",size=fontsize),
#       axis.text.x = element_text(size = 10,angle=xangle), 
#       axis.ticks.x = element_blank(),
#       axis.text.y = element_text(size=fontsize,face="plain",colour="black"),
#       legend.title= element_blank(), legend.text = element_text(size=fontsize,face="plain",colour="black"));
#   
#   return(p)
# }
# 
# 
# ## Mutect1 before
# library(ggforce)
# library(tidyr) 
# library(dplyr)
# library(ggplot2)
# library(ggforce)
# pdf(paste(OM_base,"same_order_OM_samples_data_6bases_count.pdf",sep="\\")
#     , height=12.94, width=12.94);
# xangle = 90
# fill_colors = signature_colors
# x <- mut_before_count_sum$x
# y <- mut_before_count_sum$y;
# mutation_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
# fill <- mut_before_count_sum$fill
# fill <- factor(fill, levels=mutation_types);
# samples <- unique(x);
# sample_order <- order(sapply(samples, function(s) {sum(y[which(x == s)])}), decreasing=T);
# x <- factor(x, levels=samples[sample_order]);
# data <- data.frame(x=x, y=y, fill=fill);
# 
# 
# p1 <- my_new_bar_plot(data,fill_colors=signature_colors,"Mutect before" )
# 
# 
# 
# 
# 
# ## Mutect1 after
# x <- mut_after_count_sum$x
# y <- mut_after_count_sum$y;
# mutation_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
# fill <- mut_after_count_sum$fill
# fill <- factor(fill, levels=mutation_types);
# samples <- unique(x);
# #sample_order <- order(sapply(samples, function(s) {sum(y[which(x == s)])}), decreasing=T);
# x <- factor(x, levels=samples[sample_order]);
# data <- data.frame(x=x, y=y, fill=fill);
# 
# p2 <- my_new_bar_plot(data,fill_colors=signature_colors,"Mutect After" )
# 
# grid.arrange(p1,p2,
#              nrow = 2, top = textGrob(Cancer_type,gp=gpar(fontsize=24,font=3)))
# 
# 
# dev.off()
# 
# 
# ## All the samples together (sample not in the same order)
# 
# source("C:\\Users\\abc73_000\\Documents\\GitHub\\VAF\\six_base_function_util.R")
# pdf(paste(OM_base,"OM_samples_data_6bases_count.pdf",sep="\\")
#     , height=16.94, width=12.94);
# 
# p1 <- sample_signature(mut_before_count_sum,fill_colors = signature_colors, title = "Mutect Before")
# p2 <- sample_signature(mut_after_count_sum,fill_colors = signature_colors, title = "Mutect After")
# p3 <- sample_signature(mut2_before_count_sum,fill_colors = signature_colors, title = "Mutect2 Before")
# p4 <- sample_signature(mut2_after_count_sum,fill_colors = signature_colors, title = "Mutec2 After")
# 
# 
# grid.arrange(p1,p2,
#              nrow = 2, top = textGrob(Cancer_type,gp=gpar(fontsize=24,font=3)))
# 
# grid.arrange(p3,p4,
#              nrow = 2, top = textGrob(Cancer_type,gp=gpar(fontsize=24,font=3)))
# 
# dev.off()
# 
# 
# 
# for (i in sample){
#   data <- as.data.frame(table(clean_sanger[sample_name== i,.(conver_mut_type)]))
#   colnames(data) <- c("conver_mut_type","number")
#   data <- setDT(data)
#   data_rate <- count_mutation_rates(data,signature_levels)
#   data_number <- count_mutation_number(data,signature_levels)
#   rate_data <- prepare_bar_plot_from_table(data_rate,i)
#   number_data <- prepare_bar_plot_from_table(data_number,i)
#   
#   p1 <- my_barplot(rate_data, fill_colors=signature_colors, abs_text_size=12, xangle=30)
#   p1 <- p1+labs(title="mut_Ratio", y="", x="", fill="Mut.");
#   
#   p2 <- my_barplot(number_data, fill_colors=signature_colors, abs_text_size=12, xangle=30)
#   p2 <- p2+labs(title="mut_count", y="", x="", fill="Mut.");
#   
#   
#   grid.arrange(p1,p2, ncol = 2, top = textGrob(sample,gp=gpar(fontsize=24,font=3)))
#   
# }
# dev.off()
# ###### Plot end ###
# 
# ## Analyzed the sanger Mutation vs our mutation number
# 
# 
# sanger_signature <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\OM_mutation_compare_with_Sanger\\Sanger_mutation.xlsx",
#                                sheet ="Canine", skip = 29)
# 
# clean_sanger <- sanger_signature[,c("#Chr","Position","Ref","Alt","Sample")]
# clean_sanger$chrom <- paste("chr",clean_sanger$`#Chr`,sep ="")
# clean_sanger$chrom_loc <- paste(clean_sanger$chrom,clean_sanger$Position,sep = "_")
# clean_sanger <- setDT(clean_sanger)
# 
# samples <- sort(unique(sanger_signature$Sample))
# our_data <- fread("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\OM_mutation_compare_with_Sanger\\our_totalSangerOM_mutation.txt")
# colnames(our_data) <- c("chrom","Position","Ref","Alt","Sample")
# 
# our_data_sample_convert <- sapply(our_data$Sample,convert_sample)
# our_data$Sample <- our_data_sample_convert
# our_data$chrom_loc <- paste(our_data$chrom,our_data$Position,sep ="_")
# 
# 
# our_mut_number <- c()
# sanger_mut_number <- c()
# sample_name <- c()
# 
# for (samp in samples){
#  our_num <- nrow(our_data[Sample==samp,]) 
#  sanger_num <- nrow(clean_sanger[Sample==samp,])
#  our_mut_number <- c(our_mut_number, our_num)
#  sanger_mut_number <- c(sanger_mut_number,sanger_num)
#  sample_name <- c(sample_name,samp)
#  
# }
# 
# data <- data.frame(our_mut_count= our_mut_number,
#                    sanger_mut_count= sanger_mut_number, 
#                    sample_name = sample_name )
# 
# data$diff <- abs(data$sanger_mut_count-data$our_mut_count)
# 
# highlight_df <- data %>% 
#   filter(diff >20) %>% 
#   write.table(file ="G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\OM_mutation_compare_with_Sanger\\mut_rate_most_different_sample.txt",
#               col.names = T,row.names = F,quote = F,sep ="\t")
# 
# 
# 
# 
# pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\OM_mutation_compare_with_Sanger\\Compare_OM_mut.pdf"
#     , height=4.84, width=4.84);
# 
# ggplot(data= data, aes(x=our_mut_count, y= sanger_mut_count))+
#   geom_point(shape=1,size= 3.5)+
#   geom_point(data=highlight_df, 
#              aes(x=our_mut_count,y=sanger_mut_count), 
#              color='red',
#              size=3)+
#   ylim(0,125)+
#   xlim(0,125)+
#   geom_abline(intercept = 0, slope = 1,linetype="longdash", color = "blue", size = 1)+
#   theme(axis.text=regular.text, 
#         axis.title.y=regular.text,
#         axis.title.x =regular.text,
#         axis.text.x = element_text(angle=30, hjust=1), 
#         panel.background=element_blank(), 
#         axis.line=element_line(color="black"),
#         legend.position="top", 
#         legend.title=regular.text, 
#         legend.text=regular.text, 
#         legend.key=element_blank())
# 
# dev.off()
