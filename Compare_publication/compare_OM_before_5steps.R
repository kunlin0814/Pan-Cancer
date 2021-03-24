library(tidyverse)
library(data.table)
library(gridExtra)
library(ggrepel)
library(readxl)
library(ggpubr)
library(grid)
source(#"C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
  "/Volumes/Research/GitHub/R_util/my_util.R")

source(#"C:/Users/abc73/Documents/GitHub/Breed_prediction/build_sample_meta_data.R")
  "/Volumes/Research/GitHub/Breed_prediction/build_sample_meta_data.R")
source("/Volumes/Research/GitHub/VAF/six_base_function_util.R")

whole_wes_clean_breed_table <- fread("/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt")
#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table_02_19.txt") 

fontsize=20
dot_size <- 1.4;
abs_text_size <- 16;
regular.text <- element_text(colour="black",size=abs_text_size);
xangle <- 45;
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
regular.text <- element_text(colour="black",size=20)


create_overlap_summary <- function(our_data,publish_data,intercet_sample){
  
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
    #samp = "DD0007a"
    our_each_mut <- our_data[Sample == samp, .(chrom_loc)]
    publish_each_mut <- publish_data[Sample == samp,.(chrom_loc)]
    
    
    our_each <- nrow(unique(our_data[Sample==samp,]))
    sanger_each <- nrow(unique(publish_data[Sample==samp,]))
    intercet_data <- nrow(unique(intersect(our_each_mut,publish_each_mut)))
    denominator <- our_each+sanger_each-intercet_data
    
    # count
    number_overlap <- nrow(unique(intersect(our_each_mut,publish_each_mut)))
    uniq_number_to_us <-nrow(unique(setdiff(our_each_mut,publish_each_mut))) 
    uniq_number_to_them <- nrow(unique(setdiff(publish_each_mut,our_each_mut)))
    their_mut_number <- nrow(unique(publish_each_mut))
      UGA_mut_number <- nrow(unique(our_each_mut))
    # ratio
    uniq_ratio_to_them <- nrow(unique(setdiff(publish_each_mut,our_each_mut)))/(denominator)
    uniq_ratio_to_us <- nrow(unique(setdiff(our_each_mut,publish_each_mut)))/(denominator)
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
                     uniq_num_to_publication = as.numeric(total_uniq_num_to_them),
                     uniq_num_to_uga = as.numeric(total_uniq_num_to_us),
                     share_number = as.numeric(total_share_number),
                     total_denomitor= as.numeric(total_denomitor),
                     total_their_mut_number=as.numeric(total_their_mut_number),
                     total_UGA_mut_number=as.numeric(total_UGA_mut_number))
  data <- setDT(data)
  
  return (data)
}

#### Analyzed the ratio that overlap
seperator <- "/"
base <- "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Compare_publication/OM_mutation_compare_with_Sanger"
#"G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Compare_publication\\OM_mutation_compare_with_Sanger"
original_signature <- read_excel(paste(base,"Sanger_mutation.xlsx",sep =seperator),
                                 sheet ="Canine", skip = 29)
original_signature <- setDT(original_signature)
original_signature$Chromosome <- paste("chr",original_signature$`#Chr`,sep="")
original_signature$chrom_loc <- paste(original_signature$Chromosome,original_signature$Position,sep="_")
original_signature$merge_key <- paste(original_signature$Chromosome,original_signature$Position,original_signature$Ref,original_signature$Alt,sep="_")
target_original <- original_signature[,.(Sample,Chromosome,Position,Ref,Alt,Consequence,chrom_loc,merge_key)]

#colnames(original_signature)[5] <- "sample_names"

seperator <- "/"
sanger_signature <- fread(paste(base,"DbSNP_sanger_CDS_mut_file_after_DBSNP_03_23.txt",sep =seperator))
sanger_filter_SNP <- sanger_signature

sanger_filter_SNP$chrom_loc <- paste(sanger_filter_SNP$Chromosome,sanger_filter_SNP$Position,sep = "_")
sanger_filter_SNP$merge_key <- paste(sanger_filter_SNP$Chromosome,sanger_filter_SNP$Position,sanger_filter_SNP$Ref,sanger_filter_SNP$Alt,sep="_")

samples <- sort(unique(sanger_filter_SNP$Sample))


merge_table <-  merge(x = sanger_filter_SNP, y = target_original, by.x="merge_key", by.y="merge_key", all.x = T)
#merge_table <- inner_join(sanger_filter_SNP,target_original,by= 'merge_key')
target_col <- c("Sample.x","Chromosome.x","Position.x","Ref.x","Alt.x","chrom_loc.x","Consequence")
merge_table <- merge_table[,target_col,with = F]
merge_table <- unique(merge_table)
colnames(merge_table) <- c("Sample","Chromosome","Position","Ref","Alt","chrom_loc","Consequence")
## only select SNV data ###
target_conseq <- c("missense_variant","synonymous_variant","stop_gained")

final_merge <- merge_table[Consequence %in% target_conseq,]
clean_sanger <- setDT(final_merge)

## snv file with Burair filtering ##
our_before <- fread(paste(base,"total_final_without_Gene_Burair_Filtering3_VAF_Mutect_Before_0201.txt.gz",sep = seperator));
# col_names <- c("chrom","pos","vaf","ref","alt","tRef","tAlt","sample_names",
#                "tumor_type","symbol","tumor_coverage","both_coverage")
# colnames(our_before) <- col_names
our_data <- our_before[tumor_type=="OM",.(sample_names,chrom,pos,ref,alt)]
our_data_sample_convert <- sapply(our_data$sample_names,convert_sample)
our_data$Sample <- our_data_sample_convert
our_data$chrom_loc <- paste(our_data$chrom,our_data$pos,sep ="_")
check <- our_data[Sample=="DD0066a",]

# 
# # ## indel file
# indel_file <- fread(paste(base,"total_CDS_indel_info_withGene.txt",sep =seperator))
# indel_col_names <- c("chrom","pos","ref","alt","gene_name","ensembl_id","status","sample_names")
# colnames(indel_file) <- indel_col_names
# indel_file$tumor_type <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table)
# final_indel <- indel_file[tumor_type=="OM",.(sample_names,chrom,pos,ref,alt,status)]
# final_indel$chrom_loc <- paste(final_indel$chrom,final_indel$pos,sep ="_")
# final_indel_sample_convert <- sapply(final_indel$sample_names,convert_sample)
# final_indel$Sample <- final_indel_sample_convert
# total_snv <- rbind(our_data,final_indel)

###
samples <- unique(our_data$Sample)
our_sample <- sort(unique(our_data$Sample))
their_sample <- sort(unique(clean_sanger$Sample))
intercet_sample <- intersect(their_sample,our_sample)

# pdf(paste(base,"Burair_filtering_bar_OM_Mutation_overlap_ratio.pdf",sep="\\")
#     , height=12.94, width=12.94);
png(file = paste(base,"Before_filtering_Mutation_number_compare_with_OM_publication_03_23.png",sep =seperator),
    width = 4800, height =2700, units = "px", res = 500)
data <- create_overlap_summary(our_data,clean_sanger,intercet_sample)

count_data <- melt(data, id.vars = c("sample"),
                   measure.vars= c("uniq_num_to_uga","uniq_num_to_publication","share_number"),
                   variable.name = "fill")
count_data <- count_data[order(sample)]

x <- count_data$sample
y <- count_data$value
classify <- c("uniq_num_to_uga","uniq_num_to_publication","share_number");
fill <- count_data$fill
fill <- factor(fill, levels=classify);
samples <- unique(x);
sample_order <- order(sapply(samples, function(s) {sum(y[which(x == s)])}), decreasing=T);
x <- factor(x, levels=samples[sample_order]);
plot_data <- data.frame(x=x, y=y, fill=fill);
#plot_data <- plot_data[-which(plot_data$x =="CMT-033"),]

fill_colors <- c("cyan","black","red");

p <- my_bar_function(plot_data,fill_colors = fill_colors,
                     title="UGA OM mutation number overlapped with Sanger Publication",fontsize=20)
print(p)

dev.off()

png(file = paste(base,"Before_filtering_Mutation_ratio_compare_with_OM_publication_03_23.png",sep =seperator),
    width = 4800, height =2700, units = "px", res = 500)

ratio_data <- melt(data, id.vars = c("sample"),
                   measure.vars= c("uniq_ratio_to_uga","uniq_ratio_to_publication","share_ratio"),
                   variable.name = "fill")

ratio_data <- ratio_data[order(sample)]

x <- ratio_data$sample
y <- ratio_data$value
classify <- c("uniq_ratio_to_uga","uniq_ratio_to_publication","share_ratio");
fill <- ratio_data$fill
fill <- factor(fill, levels=classify);
samples <- unique(x);
sample_order <-sample_order;
x <- factor(x, levels=samples[sample_order]);
plot_data <- data.frame(x=x, y=y, fill=fill);
#plot_data <- plot_data[-which(plot_data$x =="CMT-033"),] # remove the outlier
fill_colors <- c("cyan","black","red");

p1 <- my_bar_function(plot_data,fill_colors = fill_colors,
                      title="UGA OM mutation ratio overlapped with Sanger Publication",fontsize=20)
print(p1)
dev.off()

## Burair filtering end 


### 5steps only end
#data <- data[order(share_ratio,decreasing = T)]
#data_5setps <- data_5setps[order(share_ratio,decreasing = T)]
fwrite(data,file=paste(base,"Before_filtering_with_sanger_result.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
fwrite(our_data,file=paste(base,"before_5steps_filtering_OM_Mut_result.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
############# Plot the sanger six bases #############