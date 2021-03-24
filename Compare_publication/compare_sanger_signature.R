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
our_Burair <- fread(paste(base,"total_final_withGene_final_Filtering3_VAF_Mutect1_orientBias3_0129.gz",sep = seperator));
# col_names <- c("chrom","pos","vaf","ref","alt","tRef","tAlt","sample_names","gene_name","ensembl_id",
#                "status","tumor_type","symbol")
# colnames(our_Burair) <- col_names
our_data <- our_Burair[tumor_type=="OM",.(sample_names,chrom,pos,ref,alt,status)]

check <- our_Burair[sample_names=="CMT-33"]
a = check[vaf > 0.1,]

data <- 
hist(check$vaf,breaks = 1000)
p <- ggplot(check, aes(x=sample_names, y=vaf, color='black')) + 
  geom_jitter(size=1.6, shape=20, position=position_jitterdodge())+
  ylab("Variant allele frequency")

print(p)
dev.off()



our_data_sample_convert <- sapply(our_data$sample_names,convert_sample)
our_data$Sample <- our_data_sample_convert
our_data$chrom_loc <- paste(our_data$chrom,our_data$pos,sep ="_")


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
png(file = paste(base,"Burair_filtering_Mutation_number_compare_with_OM_publication_03_23.png",sep =seperator),
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

png(file = paste(base,"Burair_filtering_Mutation_ratio_compare_with_OM_publication_03_23.png",sep =seperator),
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

#### plot 5 steps only and no other filtering 

## snv file with 5 steps only
our_5steps <- fread(paste(base,"Total_5steps_only_without_Gene_Burair_Filtering3_VAF_Mutect__0316.txt.gz",sep = seperator));
colnames(our_5steps) <- c("chrom","pos","vaf","ref","alt","tRef","tAlt","sample_names",
                          "tumor_type","symbol")
our_5steps <- our_5steps[tumor_type=="OM",.(sample_names,chrom,pos,ref,alt)]
our_data_sample_convert <- sapply(our_5steps$sample_names,convert_sample)
our_5steps$Sample <- our_data_sample_convert
our_5steps$chrom_loc <- paste(our_5steps$chrom,our_5steps$pos,sep ="_")
samples <- unique(our_5steps$Sample)
our_sample <- sort(unique(our_5steps$Sample))
their_sample <- sort(unique(clean_sanger$Sample))
intercet_sample <- intersect(their_sample,our_sample)

# pdf(paste(base,"Burair_filtering_bar_OM_Mutation_overlap_ratio.pdf",sep="\\")
#     , height=12.94, width=12.94);
png(file = paste(base,"5steps_only_Mutation_number_compare_with_OM_publication_03_23.png",sep =seperator),
    width = 4800, height =2700, units = "px", res = 500)
data_5setps <- create_overlap_summary(our_5steps,clean_sanger,intercet_sample)

count_data <- melt(data_5setps, id.vars = c("sample"),
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

png(file = paste(base,"5steps_only_Mutation_ratio_compare_with_OM_publication_03_23.png",sep =seperator),
    width = 4800, height =2700, units = "px", res = 500)

ratio_data <- melt(data_5setps, id.vars = c("sample"),
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

### 5steps only end
#data <- data[order(share_ratio,decreasing = T)]
#data_5setps <- data_5setps[order(share_ratio,decreasing = T)]
fwrite(clean_sanger,file = paste(base,"sanger_DBSNP_CDS_mutation.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")

fwrite(our_data,file = paste(base,"Burair_filtering_Mut_results.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
fwrite(our_5steps,file=paste(base,"5steps_filtering_Mut_result.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
fwrite(data,file=paste(base,"Burair_filtering_with_sanger_result.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
fwrite(data_5setps,file = paste(base,"5steps_only_with_sanger_result.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
############# Plot the sanger six bases #############


fontsize <- 20;
signature_colors <- c("cyan","black","red","gray","green","pink");
signature_levels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
col_name <- c("number","ref","alt")
OM_base <- "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Six_base_sub_FFPE"
MT_base <- "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Six_base_sub_FFPE"
OM_sample <-  sort(list.files (path =paste(OM_base,"OM","Mutect1",sep= "\\")))
MC_sample <- sort(list.files (path = paste(MT_base,"MT","Mutect1",sep= "\\")))
fill_colors <- signature_colors 
Cancer_type <- "OM"

## to see if we want to excldued
# MC_sample <- MC_sample[-match('CMT-33', MC_sample)]
xangle <- 45

if (Cancer_type =="OM"){
  total_sample <-OM_sample
  base <- OM_base
}else{
  total_sample <- MC_sample
  base <- MT_base
}
mut_before_count_sum <- NULL
mut_after_count_sum <- NULL
mut2_before_count_sum <- NULL
mut2_after_count_sum <- NULL

sanger_signature <- read_excel(paste(base,"Sanger_mutation.xlsx",sep="\\"),
                               sheet ="Canine", skip = 29)
sanger_signature <- setDT(sanger_signature)
sample <- sort(unique(sanger_signature$Sample))
clean_sanger <- sanger_signature[, .(Ref,Alt,Sample)]
colnames(clean_sanger) <- c("ref","alt","sample_name")

clean_sanger <- clean_table(clean_sanger)

clean_sanger$conver_mut_type <- sapply(clean_sanger$mut_type, convert_mutation_type)

# pdf(paste(OM_base,"six_base_Compare_OM.pdf",sep="\\")
#     , height=8.94, width=12.84);


for (sample in total_sample){
  if(grepl("-1",sample)){
    split_words <- str_split(sample,'-')[[1]][1]
    sanger_sample <- paste(split_words,'c',sep = "")
  }else if(grepl("-2",sample)){
    split_words <- str_split(sample,'-')[[1]][1]
    sanger_sample <- paste(split_words,'d',sep = "")
  }
  else{
    sanger_sample <- paste(sample,'a',sep = "")
  }
  
  ## calculate Sanger mut and 6 bases
  # data <- as.data.frame(table(clean_sanger[sample_name== sanger_sample,.(conver_mut_type)]))
  # colnames(data) <- c("conver_mut_type","number")
  # data <- setDT(data)
  # data_rate <- count_mutation_rates(data,signature_levels)
  # data_number <- count_mutation_number(data,signature_levels)
  # rate_data <- prepare_bar_plot_from_table(data_rate,sanger_sample)
  # number_data <- prepare_bar_plot_from_table(data_number,sanger_sample)
  # 
  # ps <- my_barplot(rate_data, fill_colors=signature_colors, abs_text_size=12, xangle=30)
  # ps <- ps+labs(title="Sanger_mut_Ratio", y="", x="", fill="Mut.");
  # 
  # psn <- my_barplot(number_data, fill_colors=signature_colors, abs_text_size=12, xangle=30)
  # psn <- psn+labs(title="Sanger_mut_Count", y="", x="", fill="Mut.");
  
  ## calculate our own data
  M1before_file_name <- paste(sample,"Mutect1_before5steps.txt",sep = '_')
  M1after_file_name <- paste(sample,"Mutect1_after5steps.txt",sep = '_')
  M2before_file_name <- paste(sample,"Mutect2_before5steps.txt",sep = '_')
  M2after_file_name <- paste(sample,"Mutect2_after5steps.txt",sep = '_')
  # #
  Mutect1_before <-fread(paste(base,Cancer_type,"Mutect1",sample,"before",M1before_file_name ,sep="\\"))
  Mutect1_after <- fread(paste(base,Cancer_type,"Mutect1",sample, "after",M1after_file_name ,sep="\\"))
  Mutect2_before <-fread(paste(base,Cancer_type,"Mutect2",sample,"before",M2before_file_name ,sep="\\"))
  Mutect2_after <- fread(paste(base,Cancer_type,"Mutect2",sample, "after",M2after_file_name ,sep="\\"))
  # #
  mut_before_rate <- find_six_base(Mutect1_before,sample,rate =T);
  mut_after_rate <- find_six_base(Mutect1_after,sample,rate =T);
  mut2_before_rate <- find_six_base(Mutect2_before,sample,rate =T);
  mut2_after_rate <- find_six_base(Mutect2_after,sample,rate =T);
  # # 
  mut_before_count <- find_six_base(Mutect1_before,rate =F,sample_name=sample);
  mut_after_count <- find_six_base(Mutect1_after,rate =F,sample_name=sample);
  mut2_before_count <- find_six_base(Mutect2_before,rate =F,sample_name=sample);
  mut2_after_count <- find_six_base(Mutect2_after,rate =F,sample_name=sample);
  # # 
  mut_before_count_sum <- rbind(mut_before_count_sum,mut_before_count)
  mut_after_count_sum <- rbind(mut_after_count_sum,mut_after_count)
  mut2_before_count_sum <- rbind(mut2_before_count_sum,mut2_before_count)
  mut2_after_count_sum <- rbind(mut2_after_count_sum,mut2_after_count)
  # 
  # 
  # ## Plot the signature
    # p1 <- my_barplot(mut_before_rate, fill_colors, abs_text_size=12, xangle=30)
    # p1 <- p1 + labs(title="mut_before_Ratio", y="", x="", fill="Mut.");
    # p2 <- my_barplot(mut_after_rate,fill_colors=signature_colors, abs_text_size=12, xangle=30)
    # p2 <- p2 + labs(title="mut_after_Ratio", y="", x="", fill="Mut.");
    # p3 <- my_barplot(mut2_before_rate,fill_colors, abs_text_size=12, xangle=30)
    # p3 <- p3 + labs(title="mut2_before_Ratio", y="", x="", fill="Mut.");
    # p4 <- my_barplot(mut2_after_rate,fill_colors=signature_colors, abs_text_size=12, xangle=30)
    # p4 <- p4 + labs(title="mut2_after_Ratio", y="", x="", fill="Mut.");
  # # #   
  # # #   
    # p1_1 <- my_barplot(mut_before_count,fill_colors, abs_text_size=12, xangle=30)
    # p1_1 <- p1_1 + labs(title="mut_before_Counts", y="", x="", fill="Mut.");
    # p2_1 <- my_barplot(mut_after_count,fill_colors=signature_colors, abs_text_size=12, xangle=30)
    # p2_1 <- p2_1 + labs(title="mut_after_Counts", y="", x="", fill="Mut.");
    # p3_1 <- my_barplot(mut2_before_count,fill_colors, abs_text_size=12, xangle=30)
    # p3_1 <- p3_1 + labs(title="mut2_before_Counts", y="", x="", fill="Mut.");
    # p4_1 <- my_barplot(mut2_after_count,fill_colors=signature_colors, abs_text_size=12, xangle=30)
    # p4_1 <- p4_1 + labs(title="mut2_after_Counts", y="", x="", fill="Mut.");

  # # #   lay <- rbind(c(1,2),
  # # #                c(3,4))
  # # #   #c(6,7,8,9,9))
    # grid.arrange(p1,p2,p1_1,p2_1,
    #              p3,p4,p3_1,p4_1,
    #              nrow = 2, top = textGrob(sample,gp=gpar(fontsize=24,font=3)))

}

dev.off()

## check what samples are not normal
ex <- setDT(mut_before_count_sum)
answer <- ex[, .(total=sum(y)),keyby = .(x)][order(-total)]

### check the sample order
# mut_before_count_sum 
# mut_after_count_sum 
# mut2_before_count_sum 
# mut2_after_count_sum 

my_new_bar_plot <- function(data,fill_colors,title){
  p <- ggplot(data, aes(x=x, y=y, fill=fill)) + geom_bar(stat="identity",position='stack', width=0.6)+
    ggtitle(title)+ 
    scale_fill_manual(values=fill_colors)+
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_blank(),
      #element_text(face="plain",colour="black",size=fontsize),
      axis.title.y = element_blank(),
      #element_text(face="plain",colour="black",size=fontsize),
      axis.text.x = element_text(size = 10,angle=xangle), 
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=fontsize,face="plain",colour="black"),
      legend.title= element_blank(), legend.text = element_text(size=fontsize,face="plain",colour="black"));
  
  return(p)
}


## Mutect1 before
library(ggforce)
library(tidyr) 
library(dplyr)
library(ggplot2)
library(ggforce)
pdf(paste(OM_base,"same_order_OM_samples_data_6bases_count.pdf",sep="\\")
    , height=12.94, width=12.94);
xangle = 90
fill_colors = signature_colors
x <- mut_before_count_sum$x
y <- mut_before_count_sum$y;
mutation_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
fill <- mut_before_count_sum$fill
fill <- factor(fill, levels=mutation_types);
samples <- unique(x);
sample_order <- order(sapply(samples, function(s) {sum(y[which(x == s)])}), decreasing=T);
x <- factor(x, levels=samples[sample_order]);
data <- data.frame(x=x, y=y, fill=fill);


p1 <- my_new_bar_plot(data,fill_colors=signature_colors,"Mutect before" )





## Mutect1 after
x <- mut_after_count_sum$x
y <- mut_after_count_sum$y;
mutation_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
fill <- mut_after_count_sum$fill
fill <- factor(fill, levels=mutation_types);
samples <- unique(x);
#sample_order <- order(sapply(samples, function(s) {sum(y[which(x == s)])}), decreasing=T);
x <- factor(x, levels=samples[sample_order]);
data <- data.frame(x=x, y=y, fill=fill);

p2 <- my_new_bar_plot(data,fill_colors=signature_colors,"Mutect After" )

grid.arrange(p1,p2,
             nrow = 2, top = textGrob(Cancer_type,gp=gpar(fontsize=24,font=3)))


dev.off()


## All the samples together (sample not in the same order)

source("C:\\Users\\abc73_000\\Documents\\GitHub\\VAF\\six_base_function_util.R")
pdf(paste(OM_base,"OM_samples_data_6bases_count.pdf",sep="\\")
    , height=16.94, width=12.94);

p1 <- sample_signature(mut_before_count_sum,fill_colors = signature_colors, title = "Mutect Before")
p2 <- sample_signature(mut_after_count_sum,fill_colors = signature_colors, title = "Mutect After")
p3 <- sample_signature(mut2_before_count_sum,fill_colors = signature_colors, title = "Mutect2 Before")
p4 <- sample_signature(mut2_after_count_sum,fill_colors = signature_colors, title = "Mutec2 After")


grid.arrange(p1,p2,
             nrow = 2, top = textGrob(Cancer_type,gp=gpar(fontsize=24,font=3)))

grid.arrange(p3,p4,
             nrow = 2, top = textGrob(Cancer_type,gp=gpar(fontsize=24,font=3)))

dev.off()



for (i in sample){
  data <- as.data.frame(table(clean_sanger[sample_name== i,.(conver_mut_type)]))
  colnames(data) <- c("conver_mut_type","number")
  data <- setDT(data)
  data_rate <- count_mutation_rates(data,signature_levels)
  data_number <- count_mutation_number(data,signature_levels)
  rate_data <- prepare_bar_plot_from_table(data_rate,i)
  number_data <- prepare_bar_plot_from_table(data_number,i)
  
  p1 <- my_barplot(rate_data, fill_colors=signature_colors, abs_text_size=12, xangle=30)
  p1 <- p1+labs(title="mut_Ratio", y="", x="", fill="Mut.");
  
  p2 <- my_barplot(number_data, fill_colors=signature_colors, abs_text_size=12, xangle=30)
  p2 <- p2+labs(title="mut_count", y="", x="", fill="Mut.");
  
  
  grid.arrange(p1,p2, ncol = 2, top = textGrob(sample,gp=gpar(fontsize=24,font=3)))
  
}
dev.off()
###### Plot end ###

## Analyzed the sanger Mutation number vs our mutation number


sanger_signature <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\OM_mutation_compare_with_Sanger\\Sanger_mutation.xlsx",
                               sheet ="Canine", skip = 29)

clean_sanger <- sanger_signature[,c("#Chr","Position","Ref","Alt","Sample")]
clean_sanger$chrom <- paste("chr",clean_sanger$`#Chr`,sep ="")
clean_sanger$chrom_loc <- paste(clean_sanger$chrom,clean_sanger$Position,sep = "_")
clean_sanger <- setDT(clean_sanger)

samples <- sort(unique(sanger_signature$Sample))
our_data <- fread("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\OM_mutation_compare_with_Sanger\\our_totalSangerOM_mutation.txt")
colnames(our_data) <- c("chrom","Position","Ref","Alt","Sample")

our_data_sample_convert <- sapply(our_data$Sample,convert_sample)
our_data$Sample <- our_data_sample_convert
our_data$chrom_loc <- paste(our_data$chrom,our_data$Position,sep ="_")


our_mut_number <- c()
sanger_mut_number <- c()
sample_name <- c()

for (samp in samples){
 our_num <- nrow(our_data[Sample==samp,]) 
 sanger_num <- nrow(clean_sanger[Sample==samp,])
 our_mut_number <- c(our_mut_number, our_num)
 sanger_mut_number <- c(sanger_mut_number,sanger_num)
 sample_name <- c(sample_name,samp)
 
}

data <- data.frame(our_mut_count= our_mut_number,
                   sanger_mut_count= sanger_mut_number, 
                   sample_name = sample_name )

data$diff <- abs(data$sanger_mut_count-data$our_mut_count)

highlight_df <- data %>% 
  filter(diff >20) %>% 
  write.table(file ="G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\OM_mutation_compare_with_Sanger\\mut_rate_most_different_sample.txt",
              col.names = T,row.names = F,quote = F,sep ="\t")




pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\OM_mutation_compare_with_Sanger\\Compare_OM_mut.pdf"
    , height=4.84, width=4.84);

ggplot(data= data, aes(x=our_mut_count, y= sanger_mut_count))+
  geom_point(shape=1,size= 3.5)+
  geom_point(data=highlight_df, 
             aes(x=our_mut_count,y=sanger_mut_count), 
             color='red',
             size=3)+
  ylim(0,125)+
  xlim(0,125)+
  geom_abline(intercept = 0, slope = 1,linetype="longdash", color = "blue", size = 1)+
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =regular.text,
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.position="top", 
        legend.title=regular.text, 
        legend.text=regular.text, 
        legend.key=element_blank())

dev.off()
