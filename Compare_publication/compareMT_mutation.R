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
  #"/Volumes/Research/GitHub/VAF/six_base_function_util.R")


whole_wes_clean_breed_table <- fread(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_07.txt")
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt") 

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

# ## indel file
indel_file <- fread(paste(base,"total_CDS_indel_info_withGene_04_08.txt",sep =seperator))
indel_col_names <- c("chrom","pos","ref","alt","gene_name","ensembl_id","status","sample_names")
colnames(indel_file) <- indel_col_names
indel_file$tumor_type <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table)
final_indel <- indel_file[tumor_type=="MT",.(sample_names,chrom,pos,ref,alt,status)]
final_indel <- final_indel[-c(154:165),]
final_indel$chrom_loc <- paste(final_indel$chrom,final_indel$pos,sep ="_")
final_indel$Case <-  sapply(final_indel$sample_names,convert_MT_sample)




create_overlap_summary <- function(our_MT,publisMT,intercet_sample,method = 'min'){
  
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
    
    our_each_mut <- our_MT[Case == samp, .(chrom_loc)]
    publish_each_mut <- publisMT[Case == samp,.(chrom_loc)]
    
    
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
                     Unique_to_original_study = as.numeric(total_uniq_num_to_them),
                     Unique_to_our_study = as.numeric(total_uniq_num_to_us),
                     Shared = as.numeric(total_share_number),
                     total_denomitor= as.numeric(total_denomitor),
                     total_publication_mut_number=as.numeric(total_their_mut_number),
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


seperator = "/"
fontsize <- 20;
dot_size <- 1.4;
abs_text_size <- 16;
regular.text <- element_text(colour="black",size=abs_text_size);
xangle <- 45;
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
regular.text <- element_text(colour="black",size=20)

base <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Compare_publication/MT_mutateion_compare_with_korean"
"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Compare_publication/MT_mutateion_compare_with_korean"


##### check overlap samples ######
### analyzed Mutect2 results
#### Analyzed the ratio use bar plot before 5steps
our_MT_before <- fread(paste(base,"total_final_before_mutect2_vaf.gz",sep = seperator));
our_MT_before <- our_MT_before[symbol=="MT Korean",]

publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep=seperator))
#publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
publisMT[Case =='CMT-024']


our_MT_before$Case <- sapply(our_MT_before$sample_name, convert_MT_sample)
our_MT_before$chrom_loc <- paste(our_MT_before$chrom,our_MT_before$pos,sep= "_")
our_MT_before <- clean_table(our_MT_before)
our_MT_before$chrom_loc <- paste(our_MT_before$chrom,our_MT_before$pos,sep ="_")
fwrite(our_MT_before,file=paste(base,"04_26","MT_mutect2_before.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
our_MT_before_samples <- unique(our_MT_before$Case)
# our_sample <- unique(our_MT_before$Case)
# their_sample <- unique(publishMT$Case)
# mutect2_before_intercet_sample <- intersect(their_sample,our_sample)
# mutect2_before_diff_sample <- setdiff(their_sample,our_sample)

# 
# ### Mutect2 after 5 steps ### 
# our_MT_after <- fread(paste(base,"total_final_after_mutect2_vaf.gz",sep = seperator));
# our_MT_after <- our_MT_after[symbol=="MT Korean",]
# our_MT_after$Case <- sapply(our_MT_after$sample_name, convert_MT_sample)
# our_MT_after$chrom_loc <- paste(our_MT_after$chrom,our_MT_after$pos,sep= "_")
# our_MT_after <- clean_table(our_MT_after)
# fwrite(our_MT_before,file=paste(base,"04_26","MT_mutect2_after.txt",sep = seperator),
#        col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
# our_MT_after_samples <- unique(our_MT_after$Case)
# their_sample <- unique(publishMT$Case)
# mutect2_after_intercet_sample <- intersect(their_sample,our_sample)
# mutect2_after_diff_sample <- setdiff(their_sample,our_sample)

## Burair
our_MT_Burair <- fread(paste(base,"Final_Total_withGene_Burair_Filtering3_VAF_Mutect_orientBiasModified_04_02txt.gz",sep = seperator));
our_MT <- our_MT_Burair[tumor_type=="MT" & symbol=="MT CUK",.(sample_names,chrom,pos,ref,alt)]

our_MT$Case <- sapply(our_MT$sample_names, convert_MT_sample)
our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
fwrite(our_MT,file=paste(base,"04_26","MT_Burair_Mutect1_filtering.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
our_sample <- unique(our_MT$Case)
# their_sample <- unique(publishMT$Case)
# Burair_sample <- unique(our_MT$Case)
# Burair_diff <- setdiff(their_sample,our_sample)

total_three_intercet <- intersect(intersect(our_MT_before_samples,publishMT$Case),
                                  our_sample)

##### check overlap samples end ######
#### Main code ######
### analyzed Mutect2 results
#### Analyzed the ratio use bar plot before 5steps
our_MT_mutect2_before <- fread(paste(base,"total_final_before_mutect2_vaf.gz",sep = seperator));
our_MT_mutect2_before <- our_MT_mutect2_before[symbol=="MT Korean",]

publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep=seperator))
#publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
our_MT_mutect2_before$Case <- sapply(our_MT_mutect2_before$sample_name, convert_MT_sample)
our_MT_mutect2_before$chrom_loc <- paste(our_MT_mutect2_before$chrom,our_MT_mutect2_before$pos,sep= "_")
#our_MT_mutect2_before <- clean_table(our_MT_mutect2_before)
our_MT_mutect2_before$chrom_loc <- paste(our_MT_mutect2_before$chrom,our_MT_mutect2_before$pos,sep ="_")

png(file = paste(base,"04_26","before_UGA_mutect2_compare_count_with_MT_publication.png",sep =seperator),
    width = 5000, height =2700, units = "px", res = 500)


data <- create_overlap_summary(our_MT_mutect2_before,publisMT,total_three_intercet, method = 'min')

fwrite(data,file=paste(base,"04_26","Mutect2_before_compare_publication.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
data <- data[sample!="CMT-033"]

count_data <- melt(data, id.vars = c("sample"),
                   measure.vars= c("Unique_to_our_study","Unique_to_original_study","Shared"),
                   variable.name = "fill")
count_data <- count_data[order(sample)]

x <- count_data$sample
y <- count_data$value
classify <- c("Unique_to_our_study","Unique_to_original_study","Shared");
fill <- count_data$fill
fill <- factor(fill, levels=classify);
samples <- unique(x);
sample_order <- order(sapply(samples, function(s) {sum(y[which(x == s)])}), decreasing=T);
Mutect2_before_order <- sample_order

x <- factor(x, levels=samples[sample_order]);
plot_data <- data.frame(x=x, y=y, fill=fill);
plot_data[plot_data$x=="CMT-495",]
fill_colors <- c("cyan","black","red");
p2 <- my_bar_function(plot_data,fill_colors = fill_colors,
                      title="MuTect2 output",
                      fontsize=35)
p2 <- p2+scale_y_continuous(breaks=c(0,50,100,150))
p2 <- p2+theme(legend.position="none",
               axis.text=regular.text, 
               axis.title=regular.text)

print(p2)
dev.off()
# 
# ## Mutect2 5 step filtering
# our_MT_after <- fread(paste(base,"total_final_after_mutect2_vaf.gz",sep = seperator));
# our_MT_after <- our_MT_after[symbol=="MT Korean",]
# our_MT_after$Case <- sapply(our_MT_after$sample_name, convert_MT_sample)
# our_MT_after$chrom_loc <- paste(our_MT_after$chrom,our_MT_after$pos,sep= "_")
# our_MT_after <- clean_table(our_MT_after)
# our_MT_after$chrom_loc <- paste(our_MT_after$chrom,our_MT_after$pos,sep ="_")
# 
# 
# png(file = paste(base,"04_26","after_UGA_mutect2_compare_count_with_MT_publication.png",sep =seperator),
#     width = 5000, height =2700, units = "px", res = 500)
# 
# data <- create_overlap_summary(our_MT_after,publisMT,total_three_intercet, method = 'min')
# fwrite(data,file=paste(base,"04_26","Mutect2_after_compare_publication.txt",sep = seperator),
#        col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
# data <- data[sample!="CMT-033"]
# count_data <- melt(data, id.vars = c("sample"),
#                    measure.vars= c("Unique_to_our_study","Unique_to_original_study","Shared"),
#                    variable.name = "fill")
# count_data <- count_data[order(sample)]
# 
# x <- count_data$sample
# y <- count_data$value
# classify <- c("Unique_to_our_study","Unique_to_original_study","Shared");
# fill <- count_data$fill
# fill <- factor(fill, levels=classify);
# samples <- unique(x);
# sample_order <- Mutect2_before_order
# x <- factor(x, levels=samples[sample_order]);
# plot_data <- data.frame(x=x, y=y, fill=fill);
# fill_colors <- c("cyan","black","red");
# p2 <- my_bar_function(plot_data,fill_colors = fill_colors,
#                       title="MuTect2 output after 5-steps filtering",
#                       fontsize=35)
# p2 <- p2+scale_y_continuous(breaks=c(0,50,100,150))
# p2 <- p2+theme(legend.position="none",
#              axis.text=regular.text, 
#              axis.title=regular.text)
# #plot_data[plot_data$x=="CMT-495",]
# 
# print(p2)
# dev.off()
##### Mutect2 end #####

## Analyze MT samples with Mutect1 Burair filtering using modified
publisMT <- setDT(publishMT)
our_MT_Burair <- fread(paste(base,"Final_Total_withGene_Burair_Filtering3_VAF_Mutect_orientBiasModified_04_02txt.gz",sep = seperator));
our_MT <- our_MT_Burair[tumor_type=="MT" & symbol=="MT CUK",.(sample_names,chrom,pos,ref,alt, status)]
our_MT$Case <- sapply(our_MT$sample_names, convert_MT_sample)
our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
our_MT <- rbind(our_MT,final_indel)
png(file = paste(base,"04_26","remove_Burair_filtering_compare_with_MT_publication.png",sep =seperator),
    width = 5000, height =2700, units = "px", res = 500)

data <- create_overlap_summary(our_MT,publisMT,total_three_intercet, method ='min')
fwrite(data,file=paste(base,"04_26","Burair_filtering_compare_MT_mutect2_publication.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
data <-data[sample!="CMT-033"]
count_data <- melt(data, id.vars = c("sample"),
                   measure.vars= c("Unique_to_our_study","Unique_to_original_study","Shared"),
                   variable.name = "fill")
count_data <- count_data[order(sample)]

x <- count_data$sample
y <- count_data$value
classify <- c("Unique_to_our_study","Unique_to_original_study","Shared");
fill <- count_data$fill
fill <- factor(fill, levels=classify);
samples <- unique(x);
sample_order <-Mutect2_before_order
x <- factor(x, levels=samples[sample_order]);
plot_data <- data.frame(x=x, y=y, fill=fill);

#plot_data <- plot_data[-which(plot_data$x =="CMT-033"),]

fill_colors <- c("cyan","black","red");

p <- my_bar_function(plot_data,fill_colors = fill_colors,
                     title="Somatic mutation after 5-steps filtering & \npaired-read strand orientation filtering",fontsize=35)
p <- p+scale_y_continuous(breaks=c(0,50,100,150))
p <- p+theme(legend.position="none",
             axis.text=regular.text, 
             axis.title=regular.text)
print(p)
dev.off()

## Burair filtering end


### Analyzed with mutect1


# ## indel file
indel_file <- fread(paste(base,"total_CDS_indel_info_withGene_04_08.txt",sep =seperator))
indel_col_names <- c("chrom","pos","ref","alt","gene_name","ensembl_id","status","sample_names")
colnames(indel_file) <- indel_col_names
indel_file$tumor_type <- match_vector_table(indel_file$sample_names,"DiseaseAcronym2",whole_wes_clean_breed_table)
final_indel <- indel_file[tumor_type=="MT",.(sample_names,chrom,pos,ref,alt,status)]
final_indel <- final_indel[-c(154:165),]
final_indel$chrom_loc <- paste(final_indel$chrom,final_indel$pos,sep ="_")
final_indel$Case <-  sapply(final_indel$sample_names,convert_MT_sample)


our_MT_Mutect_before <- fread(paste(base,"total_final_before_mutect_vaf.gz",sep = seperator));
our_MT_Mutect_before <- our_MT_Mutect_before[tumor_type=="MT" & symbol=="MT Korean", .(sample_names,chrom,pos,ref,alt,status)]
# this function remove indel  but now 4/27 add strelka to replace
#our_MT_Mutect_before <- clean_table(our_MT_Mutect_before)
our_MT_Mutect_before$chrom_loc <- paste(our_MT_Mutect_before$chrom,our_MT_Mutect_before$pos,sep ="_")
our_MT_Mutect_before$Case <- sapply(our_MT_Mutect_before$sample_name, convert_MT_sample)
our_MT_Mutect_before <- rbind(our_MT_Mutect_before,final_indel)


### Mutect1 after 5 steps
our_MT_Mutect_after <- fread(paste(base,"total_final_after_mutect_vaf.gz",sep = seperator));
our_MT_Mutect_after <- our_MT_Mutect_after[tumor_type=="MT" & symbol=="MT Korean",.(sample_names,chrom,pos,ref,alt,status)]
#our_MT_Mutect_after <- clean_table(our_MT_Mutect_after)
our_MT_Mutect_after$chrom_loc <- paste(our_MT_Mutect_after$chrom,our_MT_Mutect_after$pos,sep= "_")
our_MT_Mutect_after$Case <- sapply(our_MT_Mutect_after$sample_name, convert_MT_sample)
our_MT_Mutect_after <- rbind(our_MT_Mutect_after,final_indel)


png(file = paste(base,"04_26","Mutect1_after_5steps.png",sep =seperator),
    width = 5000, height =2700, units = "px", res = 500)
data <- create_overlap_summary(our_MT_Mutect_after,publisMT,total_three_intercet, method = 'min')
fwrite(data,file=paste(base,"04_26","Mutect1_after5steps_data.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
data <- data[sample!="CMT-033"]
count_data <- melt(data, id.vars = c("sample"),
                   measure.vars= c("Unique_to_our_study","Unique_to_original_study","Shared"),
                   variable.name = "fill")
count_data <- count_data[order(sample)]

x <- count_data$sample
y <- count_data$value
classify <- c("Unique_to_our_study","Unique_to_original_study","Shared");
fill <- count_data$fill
fill <- factor(fill, levels=classify);
samples <- unique(x);
sample_order <- Mutect2_before_order
x <- factor(x, levels=samples[sample_order]);
plot_data <- data.frame(x=x, y=y, fill=fill);
fill_colors <- c("cyan","black","red");
p2 <- my_bar_function(plot_data,fill_colors = fill_colors,
                      title="Somatic mutation after 5-steps filtering",
                      fontsize=35)
p2 <- p2+scale_y_continuous(breaks=c(0,50,100,150))
p2 <- p2+theme(legend.position="none",
               axis.text=regular.text, 
               axis.title=regular.text)
#plot_data[plot_data$x=="CMT-495",]

print(p2)
dev.off()

###### Mutect1 before 5 steps ###### 
png(file = paste(base,"04_26","Mutect1_before_5steps.png",sep =seperator),
    width = 5000, height =2700, units = "px", res = 500)
data <- create_overlap_summary(our_MT_Mutect_before,publisMT,total_three_intercet, method='min')
fwrite(data,file=paste(base,"04_26","Mutect1_before5steps_data.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
data <- data[sample!="CMT-033"]
count_data <- melt(data, id.vars = c("sample"),
                   measure.vars= c("Unique_to_our_study","Unique_to_original_study","Shared"),
                   variable.name = "fill")
count_data <- count_data[order(sample)]

x <- count_data$sample
y <- count_data$value
classify <- c("Unique_to_our_study","Unique_to_original_study","Shared");
fill <- count_data$fill
fill <- factor(fill, levels=classify);
samples <- unique(x);
sample_order <- Mutect2_before_order
x <- factor(x, levels=samples[sample_order]);
plot_data <- data.frame(x=x, y=y, fill=fill);
fill_colors <- c("cyan","black","red");
p2 <- my_bar_function(plot_data,fill_colors = fill_colors,
                      title="Somatic mutation before 5-steps filtering",
                      fontsize=35)
p2 <- p2+scale_y_continuous(breaks=c(0,50,100,150))
p2 <- p2+theme(legend.position="none",
               axis.text=regular.text, 
               axis.title=regular.text)
#plot_data[plot_data$x=="CMT-495",]

print(p2)
dev.off()



### Analysed Mutect1 end ####

# ## cat all file done
# 
# 
# # after 5 steps for mutect2
# base <-   "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Compare_publication\\MT_mutateion_compare_with_korean"
# publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep="\\"))
# our_total <- after_total_summary
# #publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
# our_MT <- our_total
# publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
# publisMT <- setDT(publishMT)
# our_MT <- clean_table(our_MT)
# our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
# our_MT$sample_name <- sapply(our_MT$sample_name, convert_MT_sample)
# our_MT <- setDT(our_MT)
# samples <- unique(publishMT$Case)
# our_sample <- unique(our_MT$sample_name)
# 
# their_sample <- unique(publishMT$Case)
# intercet_sample <- intersect(their_sample,our_sample)
# 
# 
# 
# total_ratio <- NULL
# total_sample <- NULL
# for (samp in intercet_sample){
#   
#   our_each_mut <- our_MT[sample_name == samp, .(chrom_loc)]
#   mt_each_mut <- publisMT[Case == samp,.(chrom_loc)]
#   number_overlap <- nrow(intersect(our_each_mut,mt_each_mut))
#   our_each <- nrow(our_MT[sample_name==samp,])
#   mt_each <- nrow(publisMT[Case==samp,])
#   
#   overlap_ratio <- number_overlap/min(our_each,mt_each)
#   total_ratio <- c(total_ratio,overlap_ratio)
#   total_sample <- c(total_sample,samp)
#   
# }
# 
# data <- data.frame(ratio = total_ratio, sample = total_sample)
# data$tumor_type <- "MT"
# 
# pdf(paste(base,"After_filtering_MT_CDS_Mutation_overlap_ratio.pdf",sep ="\\")
#     , height=4.94, width=4.94);
# 
# p <- ggplot(data, aes(x=tumor_type, y=ratio, color='black')) + 
#   geom_jitter(size=1.6, shape=20, position=position_jitterdodge(dodge.width=3))+
#   ylab("Mutation Position Overlap Ratio")
# p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
# p <- p + scale_y_continuous(limits = c(0,1),breaks = c(0,0.5,1))
# p <- p + scale_color_manual(values = 'black')
# p <- p + theme(
#   axis.text=regular.text, 
#   axis.title=regular.text, 
#   axis.text.x = element_blank(),
#   axis.title.x =element_blank(),
#   axis.title.y=regular.text,
#   legend.position="None",
#   legend.title=regular.text, 
#   legend.text=regular.text, 
#   legend.key=element_blank(),
#   panel.background=element_blank(), 
#   axis.line=element_line(color="black"), 
#   strip.background=element_rect(color="black", fill="transparent", size=1.5), 
#   strip.text = element_text(hjust=0.5, size=12, face="plain",color="black"),)
# 
# print(p)
# 
# dev.off()
# 
# # before 5 steps for mutect2
# base <-   "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Compare_publication\\MT_mutateion_compare_with_korean"
# publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep="\\"))
# our_total <- before_total_summary
# #publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
# our_MT <-our_total
# 
# 
# publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
# publisMT <- setDT(publishMT)
# our_MT <- clean_table(our_MT)
# our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
# our_MT$sample_name <- sapply(our_MT$sample_name, convert_MT_sample)
# samples <- unique(publishMT$Case)
# our_sample <- unique(our_MT$sample_name)
# 
# their_sample <- unique(publishMT$Case)
# intercet_sample <- intersect(their_sample,our_sample)
# 
# 
# 
# total_ratio <- NULL
# total_sample <- NULL
# for (samp in intercet_sample){
#   
#   our_each_mut <- our_MT[sample_name == samp, .(chrom_loc)]
#   mt_each_mut <- publisMT[Case == samp,.(chrom_loc)]
#   number_overlap <- nrow(intersect(our_each_mut,mt_each_mut))
#   our_each <- nrow(our_MT[sample_name==samp,])
#   mt_each <- nrow(publisMT[Case==samp,])
#   
#   overlap_ratio <- number_overlap/min(our_each,mt_each)
#   total_ratio <- c(total_ratio,overlap_ratio)
#   total_sample <- c(total_sample,samp)
#   
# }
# 
# data <- data.frame(ratio = total_ratio, sample = total_sample)
# data$tumor_type <- "MT"
# 
# pdf(paste(base,"Before_filtering_MT_CDS_Mutation_overlap_ratio.pdf",sep ="\\")
#     , height=4.94, width=4.94);
# 
# p <- ggplot(data, aes(x=tumor_type, y=ratio, color='black')) + 
#   geom_jitter(size=1.6, shape=20, position=position_jitterdodge(dodge.width=3))+
#   ylab("Mutation Position Overlap Ratio")
# p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
# p <- p + scale_y_continuous(limits = c(0,1),breaks = c(0,0.5,1))
# p <- p + scale_color_manual(values = 'black')
# p <- p + theme(
#   axis.text=regular.text, 
#   axis.title=regular.text, 
#   axis.text.x = element_blank(),
#   axis.title.x =element_blank(),
#   axis.title.y=regular.text,
#   legend.position="None",
#   legend.title=regular.text, 
#   legend.text=regular.text, 
#   legend.key=element_blank(),
#   panel.background=element_blank(), 
#   axis.line=element_line(color="black"), 
#   strip.background=element_rect(color="black", fill="transparent", size=1.5), 
#   strip.text = element_text(hjust=0.5, size=12, face="plain",color="black"),)
# 
# print(p)
# 
# dev.off()
# 
# ## Analyzed used Mutect1
# source("C:\\Users\\abc73_000\\Documents\\GitHub\\VAF\\six_base_function_util.R")
# 
# excldue <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan-Cancer-Manuscript\\Figure1\\Original_Data_summary.xlsx",
#                       sheet = "Before_Matching_excluded")
# 
# create_overlap_summary <- function(our_MT,publisMT,intercet_sample){
#   
#   total_uniq_num_to_them <- NULL
#   total_uniq_num_to_us <- NULL
#   total_share_number <- NULL
#   total_uniq_ratio_to_them <- NULL
#   total_uniq_ratio_to_us <- NULL
#   total_share_ratio <- NULL
#   total_sample <- NULL
#   total_denomitor <- NULL
#   for (samp in intercet_sample){
#     
#     our_each_mut <- our_MT[sample_name == samp, .(chrom_loc)]
#     publish_each_mut <- publisMT[Case == samp,.(chrom_loc)]
#     
#     
#     our_each <- nrow(our_MT[sample_name==samp,])
#     sanger_each <- nrow(publisMT[Case==samp,])
#     intercet_data <- nrow(intersect(our_each_mut,publish_each_mut))
#     denominator <- our_each+sanger_each-intercet_data
#     
#     # count
#     number_overlap <- nrow(intersect(our_each_mut,publish_each_mut))
#     uniq_number_to_us <-nrow(setdiff(our_each_mut,publish_each_mut)) 
#     uniq_number_to_them <- nrow(setdiff(publish_each_mut,our_each_mut)) 
#     
#     # ratio
#     uniq_ratio_to_them <- nrow(setdiff(publish_each_mut,our_each_mut))/(denominator)
#     uniq_ratio_to_us <- nrow(setdiff(our_each_mut,publish_each_mut))/(denominator)
#     overlap_ratio <- number_overlap/(denominator)
#     
#     
#     # summary
#     total_share_ratio <- c(total_share_ratio,overlap_ratio)
#     total_uniq_ratio_to_us <- c(total_uniq_ratio_to_us,uniq_ratio_to_us)
#     total_uniq_ratio_to_them <- c(total_uniq_ratio_to_them,uniq_ratio_to_them)
#     total_uniq_num_to_them <- c(total_uniq_num_to_them,uniq_number_to_them)
#     total_uniq_num_to_us <- c(total_uniq_num_to_us,uniq_number_to_us)
#     total_share_number <- c(total_share_number,number_overlap)
#     total_sample <- c(total_sample,samp)
#     total_denomitor <- c(total_denomitor,denominator)
#     
#   }
#   
#   data <- data.frame(share_ratio = as.numeric(total_share_ratio), 
#                      sample = total_sample,
#                      uniq_ratio_to_uga = as.numeric(total_uniq_ratio_to_us),
#                      uniq_ratio_to_publication = as.numeric(total_uniq_ratio_to_them),
#                      uniq_num_to_publication = as.numeric(total_uniq_num_to_them),
#                      uniq_num_to_uga = as.numeric(total_uniq_num_to_us),
#                      share_number = as.numeric(total_share_number),
#                      total_denomitor= as.numeric(total_denomitor))
#   data <- setDT(data)
#   
#   return (data)
# }
# 
# 
# 
# 
# 
# ## Anayzed used Mutect2 before mutation Distribution
# base <- "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Compare_publication\\MT_mutateion_compare_with_korean"
# 
# publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep="\\"))
# our_total <- before_total_summary
# 
# #publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
# our_MT <- our_total
# publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
# publisMT <- setDT(publishMT)
# our_MT <- clean_table(our_MT)
# our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
# our_MT$sample_name <- sapply(our_MT$sample_name, convert_MT_sample)
# samples <- unique(publishMT$Case)
# 
# 
# our_mut_number <- c()
# pub_mut_number <- c()
# sample_name <- c()
# 
# for (samp in samples){
#   our_num <- nrow(our_total[sample_name==samp,]) 
#   pb_num <- nrow(publishMT[Case==samp,])
#   our_mut_number <- c(our_mut_number, our_num)
#   pub_mut_number <- c(pub_mut_number,pb_num)
#   sample_name <- c(sample_name,samp)
#   
# }
# 
# data <- data.frame(our_mut_count= our_mut_number,
#                    pub_mut_number= pub_mut_number, 
#                    sample_name = sample_name )
# 
# 
# 
# pdf(paste(base,"Before_filtering_Compare_MT_mut_mutect2.pdf",sep="\\")
#     , height=4.84, width=4.84);
# 
# ggplot(data= data, aes(x=our_mut_count, y= pub_mut_number))+
#   geom_point(shape=1,size= 3.5)+
#   # geom_point(data=highlight_df, 
#   #            aes(x=our_mut_count,y=sanger_mut_count), 
#   #            color='red',
#   #            size=3)+
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
# 
# 
# ## Anayzed used Mutect2 after mutation distribution
# 
# publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep="\\"))
# our_total <- after_total_summary
# 
# #publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
# our_MT <- our_total
# publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
# publisMT <- setDT(publishMT)
# our_MT <- clean_table(our_MT)
# our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
# our_MT$sample_name <- sapply(our_MT$sample_name, convert_MT_sample)
# samples <- unique(publishMT$Case)
# 
# 
# our_mut_number <- c()
# pub_mut_number <- c()
# sample_name <- c()
# 
# for (samp in samples){
#   our_num <- nrow(our_total[sample_name==samp,]) 
#   pb_num <- nrow(publishMT[Case==samp,])
#   our_mut_number <- c(our_mut_number, our_num)
#   pub_mut_number <- c(pub_mut_number,pb_num)
#   sample_name <- c(sample_name,samp)
#   
# }
# 
# data <- data.frame(our_mut_count= our_mut_number,
#                    pub_mut_number= pub_mut_number, 
#                    sample_name = sample_name )
# 
# 
# 
# pdf(paste(base,"After_filtering_Compare_MT_mut_mutect2.pdf",sep="\\")
#     , height=4.84, width=4.84);
# 
# ggplot(data= data, aes(x=our_mut_count, y= pub_mut_number))+
#   geom_point(shape=1,size= 3.5)+
#   # geom_point(data=highlight_df, 
#   #            aes(x=our_mut_count,y=sanger_mut_count), 
#   #            color='red',
#   #            size=3)+
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
# 
# 
# 
# ## Analyzed Mut1 distribution
# 
# ## after
# publishMT <- fread(paste(base,"publish_MT_cds.txt",sep="\\"))
# our_total <- fread(paste(base,"After_Total_VAF_depth_summary.txt",sep = "\\"));
# 
# publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
# our_MT <- our_total[tumor_type=="MT",]
# 
# 
# publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
# publisMT <- setDT(publishMT)
# our_MT <- clean_table(our_MT)
# our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
# our_MT$Case <- sapply(our_MT$sample_name, convert_MT_sample)
# samples <- unique(publishMT$Case)
# 
# 
# our_mut_number <- c()
# pub_mut_number <- c()
# sample_name <- c()
# 
# for (samp in samples){
#   our_num <- nrow(our_total[sample_name==samp,]) 
#   pb_num <- nrow(publishMT[Case==samp,])
#   our_mut_number <- c(our_mut_number, our_num)
#   pub_mut_number <- c(pub_mut_number,pb_num)
#   sample_name <- c(sample_name,samp)
#   
# }
# 
# data <- data.frame(our_mut_count= our_mut_number,
#                    pub_mut_number= pub_mut_number, 
#                    sample_name = sample_name )
# 
# 
# 
# pdf(paste(base,"After_Compare_MT_mut.pdf",sep="\\")
#     , height=4.84, width=4.84);
# 
# ggplot(data= data, aes(x=our_mut_count, y= pub_mut_number))+
#   geom_point(shape=1,size= 3.5)+
#   # geom_point(data=highlight_df, 
#   #            aes(x=our_mut_count,y=sanger_mut_count), 
#   #            color='red',
#   #            size=3)+
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
# ## before mutect1
# publishMT <- fread(paste(base,"publish_MT_cds.txt",sep="\\"))
# our_total <- fread(paste(base,"Before_Total_VAF_sum.txt",sep = "\\"));
# 
# publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
# our_MT <- our_total[tumor_type=="MT",]
# 
# 
# publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
# publisMT <- setDT(publishMT)
# our_MT <- clean_table(our_MT)
# our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
# our_MT$Case <- sapply(our_MT$sample_name, convert_MT_sample)
# samples <- unique(publishMT$Case)
# 
# 
# our_mut_number <- c()
# pub_mut_number <- c()
# sample_name <- c()
# 
# for (samp in samples){
#   our_num <- nrow(our_total[sample_name==samp,]) 
#   pb_num <- nrow(publishMT[Case==samp,])
#   our_mut_number <- c(our_mut_number, our_num)
#   pub_mut_number <- c(pub_mut_number,pb_num)
#   sample_name <- c(sample_name,samp)
#   
# }
# 
# data <- data.frame(our_mut_count= our_mut_number,
#                    pub_mut_number= pub_mut_number, 
#                    sample_name = sample_name )
# 
# 
# 
# pdf(paste(base,"Before_Compare_MT_mut.pdf",sep="\\")
#     , height=4.84, width=4.84);
# 
# ggplot(data= data, aes(x=our_mut_count, y= pub_mut_number))+
#   geom_point(shape=1,size= 3.5)+
#   # geom_point(data=highlight_df, 
#   #            aes(x=our_mut_count,y=sanger_mut_count), 
#   #            color='red',
#   #            size=3)+
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
# 
# 
# ## Analyzed MT six bases before and after
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
# source("C:\\Users\\abc73_000\\Documents\\GitHub\\Pan-Cancer\\six_base_function_util.R")
# 
# regular.text <- element_text(colour="black",size=20)
# fontsize <- 20;
# xangle <- 90
# signature_colors <- c("cyan","black","red","gray","green","pink");
# signature_levels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
# base <- "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\VAF"
# 
# before_total <- fread(paste(base,"Before_Total_VAF_sum.txt",sep="\\"))
# after_total <- fread(paste(base,"After_Total_VAF_sum.txt",sep="\\"))
# 
# before_MT <- before_total[tumor_type=="MT",]
# after_MT <- after_total[tumor_type=="MT",]
# 
# total_sample <- sort(unique(before_MT$sample_name))
# 
# mut_before_count_sum <- NULL
# mut_after_count_sum <- NULL
# 
# for (sample in total_sample){
#   before <- conver_to_six_bases_basedon_mutation(before_MT,sample,"MT",rate=F)
#   after <- conver_to_six_bases_basedon_mutation(after_MT,sample,"MT",rate = F)
#   
#   mut_before_count_sum <- rbind(mut_before_count_sum,before)
#   mut_after_count_sum <- rbind(mut_after_count_sum,after)
# }
# 
# 
# 
# mut_before_count_sum <- mut_before_count_sum[x!="CMT-33",]
# mut_after_count_sum <- mut_after_count_sum[x!="CMT-33",]
# ## Mutect1 before
# 
# library(ggforce)
# library(tidyr) 
# library(dplyr)
# library(ggplot2)
# library(ggforce)
# 
# pdf(paste(base,"MT_samples_data_6bases_count.pdf",sep="\\")
#     , height=16.94, width=12.94);
# 
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
#              nrow = 2)
# 
# 
# dev.off()
# 
# 
# 
# # ### cat all of Mutect2 data together because the publication used Mutect2
# 
# 
# cancer_type <- "MT"
# base_dir = #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate/VAF/Mutect1"
#   "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Compare_publication/MT_mutateion_compare_with_korean/Mutect2/MT"
# 
# # # coverage_file <- read_excel(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate/New_WES_QC_dataset.xlsx",
# # #   "C:\\Users\\abc73_000\\Desktop\\New_WES_QC_dataset.xlsx",
# # #   sheet = "Mut")
# #
# # coverage_file <- setDT(coverage_file)
# 
# before_total_summary <- NULL
# after_total_summary <- NULL
# 
# for (cancer in cancer_type){
# 
#   sample_names <- list.files(base_dir)
# 
#   for (sample in sample_names){
# 
#     M2before_file_name <- paste(sample,"vaf_before.txt",sep = '_')
#     M2after_file_name <- paste(sample,"vaf_after.txt",sep = '_')
# 
#     Mutect2_before <-fread(paste(base_dir,sample,M2before_file_name ,sep=seperator))
#     Mutect2_after <- fread(paste(base_dir,sample,M2after_file_name ,sep=seperator))
#     # before
#     if (nrow(Mutect2_before)!=0){
#       before_each_file <-fread(paste(base_dir,sample,M2before_file_name ,sep=seperator))
#       colnames(before_each_file) <- c("chrom","pos","vaf","ref","alt","tRef","tAlt","sample_name")
#       before_each_file$sample_name <-sample
#     }
# 
# 
#     if (nrow(Mutect2_after)!=0){
#       
#       # after
#       after_each_file <- fread(paste(base_dir,sample,M2after_file_name ,sep=seperator))
#       colnames(after_each_file) <- c("chrom","pos","vaf","ref","alt","tRef","tAlt","sample_name")
#       after_each_file$sample_name <-sample
# 
#     }
#     before_total_summary <- rbind(before_total_summary,before_each_file)
#     after_total_summary <- rbind(after_total_summary,after_each_file)
#   }
# }
# 
# before_total_summary$sample_name
# ## cat all file and write output
# 
# fwrite(before_total_summary,file = paste(base_dir,"MT_mutect2_before.txt",sep = seperator),
#        sep="\t", col.names = T, row.names = F, quote = F)
# 
# 
# fwrite(after_total_summary,file = paste(base_dir,"MT_mutect2_after.txt",sep = seperator),
#        sep="\t", col.names = T, row.names = F, quote = F)
# 
# 
# unique(before_total_summary$sample_name)
