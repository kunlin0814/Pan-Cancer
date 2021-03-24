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

exclude <- unique(unlist(whole_wes_clean_breed_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

create_overlap_summary <- function(our_MT,publisMT,intercet_sample){
  
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
    
    
    our_each <- nrow(our_MT[Case==samp,])
    sanger_each <- nrow(publisMT[Case==samp,])
    intercet_data <- nrow(intersect(our_each_mut,publish_each_mut))
    denominator <- our_each+sanger_each-intercet_data
    
    # count
    number_overlap <- nrow(intersect(our_each_mut,publish_each_mut))
    uniq_number_to_us <-nrow(setdiff(our_each_mut,publish_each_mut)) 
    uniq_number_to_them <- nrow(setdiff(publish_each_mut,our_each_mut)) 
    their_mut_number <- nrow(unique(publish_each_mut))
    UGA_mut_number <- nrow(unique(our_each_mut))
    # ratio
    uniq_ratio_to_them <- nrow(setdiff(publish_each_mut,our_each_mut))/(denominator)
    uniq_ratio_to_us <- nrow(setdiff(our_each_mut,publish_each_mut))/(denominator)
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

seperator = "/"
fontsize <- 20;
dot_size <- 1.4;
abs_text_size <- 16;
regular.text <- element_text(colour="black",size=abs_text_size);
xangle <- 45;
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
regular.text <- element_text(colour="black",size=20)

## Analyze MT samples with mutect2 and Burair filtering

base <- "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Compare_publication/MT_mutateion_compare_with_korean"
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Compare_publication/MT_mutateion_compare_with_korean"
publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep=seperator))

our_MT_Burair <- fread(paste(base,"total_final_withGene_final_Filtering3_VAF_Mutect1_orientBias3_0129.gz",sep = seperator));
our_MT <- our_MT_Burair[tumor_type=="MT" & symbol=="MT CUK",.(sample_names,chrom,pos,ref,alt)]
#colnames(our_MT) <- colnames(publishMT)
# 
our_MT$sample_names
our_MT_after <- fread(paste(base,"MT_mutect2_after.txt",sep = seperator));
our_MT_before <- fread(paste(base,"MT_mutect2_before.txt",sep = seperator));
#publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
our_MT$Case <- sapply(our_MT$sample_names, convert_MT_sample)
our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")

our_MT_after <- clean_table(our_MT_after)
our_MT_after$chrom_loc <- paste(our_MT_after$chrom,our_MT_after$pos,sep= "_")
our_MT_after$sample_name <- sapply(our_MT_after$sample_name, convert_MT_sample)
our_MT_before <- clean_table(our_MT_before)
our_MT_before$chrom_loc <- paste(our_MT_before$chrom,our_MT_before$pos,sep ="_")
our_MT_before$sample_name <- sapply(our_MT_before$sample_name, convert_MT_sample)

# colnames(our_MT_before) 
# colnames(our_MT_after) <- colnames(publishMT)

samples <- unique(our_MT$Case)
our_sample <- unique(our_MT$Case)
their_sample <- unique(publishMT$Case)
intercet_sample <- intersect(their_sample,our_sample)
setdiff(their_sample,our_sample)


## Analyzed ratio use bar plot

#### Analyzed the ratio use bar plot after 5steps
# pdf(paste(base,"Burair_filtering_bar_MT_Mutation_overlap_ratio_for_mutect2.pdf",sep="\\")
#     , height=12.94, width=12.94);

png(file = paste(base,"Burair_filtering_compare_with_MT_publication.png",sep =seperator),
    width = 4800, height =2700, units = "px", res = 500)

data <- create_overlap_summary(our_MT,publisMT,intercet_sample)

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
                     title="UGA MT mutation number overlapped with MT CUK mutect2",fontsize=20)
print(p)
dev.off()

png(file = paste(base,"Burair_filtering_compare_with_MT_publication_ratio.png",sep =seperator),
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
plot_data <- plot_data[-which(plot_data$x =="CMT-033"),] # remove the outlier
fill_colors <- c("cyan","black","red");

p1 <- my_bar_function(plot_data,fill_colors = fill_colors,
                     title="UGA MT mutation ratio overlapped with MT CUK mutect2",fontsize=20)
print(p1)
dev.off()

fwrite(our_MT,file=paste(base,"Our_MT_mut_Bur_filtering.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")
fwrite(data,file=paste(base,"Burair_filtering_with_MT_CDS_publication_result.txt",sep = seperator),
       col.names = T,row.names = F,quote = F, eol = "\n",sep = "\t")

## Burair filtering end

### analyzed Mutect2 results
#### Analyzed the ratio use bar plot before 5steps
our_MT_after <- fread(paste(base,"MT_mutect2_after.txt",sep = seperator));
our_MT_before <- fread(paste(base,"MT_mutect2_before.txt",sep = seperator));
#publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
our_MT_before$Case <- sapply(our_MT_before$sample_name, convert_MT_sample)
our_MT_before$chrom_loc <- paste(our_MT_before$chrom,our_MT_before$pos,sep= "_")
our_MT_before <- clean_table(our_MT_before)
our_MT_before$chrom_loc <- paste(our_MT_before$chrom,our_MT_before$pos,sep ="_")


# colnames(our_MT_before) 
# colnames(our_MT_after) <- colnames(publishMT)

samples <- unique(our_MT_before$Case)
our_sample <- unique(our_MT_before$Case)
their_sample <- unique(publishMT$Case)
intercet_sample <- intersect(their_sample,our_sample)
setdiff(their_sample,our_sample)


png(file = paste(base,"UGA_mutect2_compare_count_with_MT_publication.png",sep =seperator),
    width = 4800, height =2700, units = "px", res = 500)


data <- create_overlap_summary(our_MT_before,publisMT,intercet_sample)

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
fill_colors <- c("cyan","black","red");
p2 <- my_bar_function(plot_data,fill_colors = fill_colors,
                     title="Before MT mutation number overlapped with MT Korean mutect2",
                     fontsize=20)
print(p2)
dev.off()
png(file = paste(base,"UGA_mutect2_compare_ratio_with_MT_publication.png",sep =seperator),
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
fill_colors <- c("cyan","black","red");

p3 <- my_bar_function(plot_data,fill_colors = fill_colors,
                     title="Before MT mutation ratio overlapped with MT Korean mutect2",
                     fontsize=20)
print(p3)
dev.off()
##### Mutect2 end #####


### Analyzed with mutect1
our_MT_after <- fread(paste(base,"total_final_after_mutect_vaf.gz",sep = seperator));
our_MT_before <- fread(paste(base,"total_final_before_mutect_vaf.gz",sep = seperator));

our_MT_before <- our_MT_before[tumor_type=="MT" & symbol=="MT Korean", .(sample_names,chrom,pos,ref,alt)]
our_MT_after <- our_MT_after[tumor_type=="MT" & symbol=="MT Korean",.(sample_names,chrom,pos,ref,alt)]
publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
our_MT_after <- clean_table(our_MT_after)
our_MT_after$chrom_loc <- paste(our_MT_after$chrom,our_MT_after$pos,sep= "_")
our_MT_after$sample_name <- sapply(our_MT_after$sample_name, convert_MT_sample)
our_MT_before <- clean_table(our_MT_before)
our_MT_before$chrom_loc <- paste(our_MT_before$chrom,our_MT_before$pos,sep ="_")
our_MT_before$sample_name <- sapply(our_MT_before$sample_name, convert_MT_sample)

samples <- unique(publishMT$Case)
our_sample <- unique(our_MT_after$sample_name)
their_sample <- unique(publishMT$Case)
intercet_sample <- intersect(their_sample,our_sample)
setdiff(their_sample,our_sample)


data <- create_overlap_summary(our_MT_after,publisMT,intercet_sample)

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
fill_colors <- c("cyan","black","red");

p4 <- my_bar_function(plot_data,fill_colors = fill_colors,
                     title="After MT mutation number overlapped with MT Korean mutect1",fontsize=20)
print(p4)

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
fill_colors <- c("cyan","black","red");

p5 <- my_bar_function(plot_data,fill_colors = fill_colors,
                      title="After MT mutation ratio overlapped with MT Korean mutect1",fontsize=20)
print(p5)

#### Analyzed the ratio use bar plot before 5steps

data <- create_overlap_summary(our_MT_before,publisMT,intercet_sample)

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
fill_colors <- c("cyan","black","red");
p6 <- my_bar_function(plot_data,fill_colors = fill_colors,
                      title="Before MT mutation number overlapped with MT Korean mutect1",
                      fontsize=20)
print(p6)

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
fill_colors <- c("cyan","black","red");

p7 <- my_bar_function(plot_data,fill_colors = fill_colors,
                      title="Before MT mutation ratio overlapped with MT Korean mutect1",
                      fontsize=20)
print(p7)

#grid.arrange(p, p2,p3,p4,pnrow = 2)
#grid.arrange(p, p2,p3,p4,pnrow = 2)
dev.off()

### Bar plot end ##

## Anayzed ratio old ###
total_ratio <- NULL
total_sample <- NULL
for (samp in intercet_sample){
  
  our_each_mut <- our_MT[sample_name == samp, .(chrom_loc)]
  mt_each_mut <- publisMT[Case == samp,.(chrom_loc)]
  number_overlap <- nrow(intersect(our_each_mut,mt_each_mut))
  our_each <- nrow(our_MT[sample_name==samp,])
  mt_each <- nrow(publisMT[Case==samp,])
  
  overlap_ratio <- number_overlap/min(our_each,mt_each)
  total_ratio <- c(total_ratio,overlap_ratio)
  total_sample <- c(total_sample,samp)
  
}

data <- data.frame(ratio = total_ratio, sample = total_sample)
data$tumor_type <- "MT"

pdf(paste(base,"After_filtering_MT_CDS_Mutation_overlap_ratio.pdf",sep ="\\")
    , height=4.94, width=4.94);

p <- ggplot(data, aes(x=tumor_type, y=ratio, color='black')) + 
  geom_jitter(size=1.6, shape=20, position=position_jitterdodge(dodge.width=3))+
  ylab("Mutation Position Overlap Ratio")
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_y_continuous(limits = c(0,1),breaks = c(0,0.5,1))
p <- p + scale_color_manual(values = 'black')
p <- p + theme(
  axis.text=regular.text, 
  axis.title=regular.text, 
  axis.text.x = element_blank(),
  axis.title.x =element_blank(),
  axis.title.y=regular.text,
  legend.position="None",
  legend.title=regular.text, 
  legend.text=regular.text, 
  legend.key=element_blank(),
  panel.background=element_blank(), 
  axis.line=element_line(color="black"), 
  strip.background=element_rect(color="black", fill="transparent", size=1.5), 
  strip.text = element_text(hjust=0.5, size=12, face="plain",color="black"),)

print(p)

dev.off()


### cat all of Mutect2 data together

cancer_type <- "MT"
base_dir = #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate/VAF/Mutect1"
  "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Compare_publication\\MT_mutateion_compare_with_korean\\Mutect2\\MT"

# # coverage_file <- read_excel(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate/New_WES_QC_dataset.xlsx",
# #   "C:\\Users\\abc73_000\\Desktop\\New_WES_QC_dataset.xlsx",
# #   sheet = "Mut")
# 
# coverage_file <- setDT(coverage_file)

before_total_summary <- NULL
after_total_summary <- NULL

for (cancer in cancer_type){
  
  sample_names <- list.files(base_dir)
  
  for (sample in sample_names){
    
    M2before_file_name <- paste(sample,"vaf_before.txt",sep = '_')
    M2after_file_name <- paste(sample,"vaf_after.txt",sep = '_')
    
    Mutect2_before <-fread(paste(base_dir,sample,M2before_file_name ,sep="\\"))
    Mutect2_after <- fread(paste(base_dir,sample,M2after_file_name ,sep="\\"))
    
    
    
    if (nrow(Mutect2_after)!=0){
      #print(file_loc)
      before_each_file <-fread(paste(base_dir,sample,M2before_file_name ,sep="\\"))
      after_each_file <- fread(paste(base_dir,sample,M2after_file_name ,sep="\\"))
      # before
      colnames(before_each_file) <- c("chrom","pos","vaf","ref","alt")
      before_each_file$sample_name <-sample
      #coverage_info <- coverage_file[Case_ID==sample,.(coverage_stat,tumor_coverage,both_status)]
      # before_each_file$tumor_coverage <- coverage_info$tumor_coverage[1]
      # before_each_file$both_coverage <- coverage_info$both_status[1]
      # before_each_file$sample_mean <- coverage_info$coverage_stat[1]
      # before_each_file$tumor_type <- cancer
      
      # after
      colnames(after_each_file) <- c("chrom","pos","vaf","ref","alt")
      after_each_file$sample_name <-sample
      
      #coverage_info <- coverage_file[Case_ID==sample,.(coverage_stat,tumor_coverage,both_status)]
      # after_each_file$tumor_coverage <- coverage_info$tumor_coverage[1]
      # after_each_file$both_coverage <- coverage_info$both_status[1]
      # after_each_file$sample_mean <- coverage_info$coverage_stat[1]
      # after_each_file$tumor_type <- cancer
      
    }
    before_total_summary <- rbind(before_total_summary,before_each_file)
    after_total_summary <- rbind(after_total_summary,after_each_file)
  }
  #total_summary <- rbind(total_summary,after_each_file)
}
## cat all file done


# after 5 steps for mutect2
base <-   "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Compare_publication\\MT_mutateion_compare_with_korean"
publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep="\\"))
our_total <- after_total_summary
#publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
our_MT <- our_total
publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
our_MT <- clean_table(our_MT)
our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
our_MT$sample_name <- sapply(our_MT$sample_name, convert_MT_sample)
our_MT <- setDT(our_MT)
samples <- unique(publishMT$Case)
our_sample <- unique(our_MT$sample_name)

their_sample <- unique(publishMT$Case)
intercet_sample <- intersect(their_sample,our_sample)



total_ratio <- NULL
total_sample <- NULL
for (samp in intercet_sample){
  
  our_each_mut <- our_MT[sample_name == samp, .(chrom_loc)]
  mt_each_mut <- publisMT[Case == samp,.(chrom_loc)]
  number_overlap <- nrow(intersect(our_each_mut,mt_each_mut))
  our_each <- nrow(our_MT[sample_name==samp,])
  mt_each <- nrow(publisMT[Case==samp,])
  
  overlap_ratio <- number_overlap/min(our_each,mt_each)
  total_ratio <- c(total_ratio,overlap_ratio)
  total_sample <- c(total_sample,samp)
  
}

data <- data.frame(ratio = total_ratio, sample = total_sample)
data$tumor_type <- "MT"

pdf(paste(base,"After_filtering_MT_CDS_Mutation_overlap_ratio.pdf",sep ="\\")
    , height=4.94, width=4.94);

p <- ggplot(data, aes(x=tumor_type, y=ratio, color='black')) + 
  geom_jitter(size=1.6, shape=20, position=position_jitterdodge(dodge.width=3))+
  ylab("Mutation Position Overlap Ratio")
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_y_continuous(limits = c(0,1),breaks = c(0,0.5,1))
p <- p + scale_color_manual(values = 'black')
p <- p + theme(
  axis.text=regular.text, 
  axis.title=regular.text, 
  axis.text.x = element_blank(),
  axis.title.x =element_blank(),
  axis.title.y=regular.text,
  legend.position="None",
  legend.title=regular.text, 
  legend.text=regular.text, 
  legend.key=element_blank(),
  panel.background=element_blank(), 
  axis.line=element_line(color="black"), 
  strip.background=element_rect(color="black", fill="transparent", size=1.5), 
  strip.text = element_text(hjust=0.5, size=12, face="plain",color="black"),)

print(p)

dev.off()

# before 5 steps for mutect2
base <-   "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Compare_publication\\MT_mutateion_compare_with_korean"
publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep="\\"))
our_total <- before_total_summary
#publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
our_MT <-our_total


publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
our_MT <- clean_table(our_MT)
our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
our_MT$sample_name <- sapply(our_MT$sample_name, convert_MT_sample)
samples <- unique(publishMT$Case)
our_sample <- unique(our_MT$sample_name)

their_sample <- unique(publishMT$Case)
intercet_sample <- intersect(their_sample,our_sample)



total_ratio <- NULL
total_sample <- NULL
for (samp in intercet_sample){
  
  our_each_mut <- our_MT[sample_name == samp, .(chrom_loc)]
  mt_each_mut <- publisMT[Case == samp,.(chrom_loc)]
  number_overlap <- nrow(intersect(our_each_mut,mt_each_mut))
  our_each <- nrow(our_MT[sample_name==samp,])
  mt_each <- nrow(publisMT[Case==samp,])
  
  overlap_ratio <- number_overlap/min(our_each,mt_each)
  total_ratio <- c(total_ratio,overlap_ratio)
  total_sample <- c(total_sample,samp)
  
}

data <- data.frame(ratio = total_ratio, sample = total_sample)
data$tumor_type <- "MT"

pdf(paste(base,"Before_filtering_MT_CDS_Mutation_overlap_ratio.pdf",sep ="\\")
    , height=4.94, width=4.94);

p <- ggplot(data, aes(x=tumor_type, y=ratio, color='black')) + 
  geom_jitter(size=1.6, shape=20, position=position_jitterdodge(dodge.width=3))+
  ylab("Mutation Position Overlap Ratio")
p <- p + stat_summary(fun=median, fun.min=median, fun.max=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_y_continuous(limits = c(0,1),breaks = c(0,0.5,1))
p <- p + scale_color_manual(values = 'black')
p <- p + theme(
  axis.text=regular.text, 
  axis.title=regular.text, 
  axis.text.x = element_blank(),
  axis.title.x =element_blank(),
  axis.title.y=regular.text,
  legend.position="None",
  legend.title=regular.text, 
  legend.text=regular.text, 
  legend.key=element_blank(),
  panel.background=element_blank(), 
  axis.line=element_line(color="black"), 
  strip.background=element_rect(color="black", fill="transparent", size=1.5), 
  strip.text = element_text(hjust=0.5, size=12, face="plain",color="black"),)

print(p)

dev.off()

## Analyzed used Mutect1
source("C:\\Users\\abc73_000\\Documents\\GitHub\\VAF\\six_base_function_util.R")

excldue <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan-Cancer-Manuscript\\Figure1\\Original_Data_summary.xlsx",
                      sheet = "Before_Matching_excluded")

create_overlap_summary <- function(our_MT,publisMT,intercet_sample){
  
  total_uniq_num_to_them <- NULL
  total_uniq_num_to_us <- NULL
  total_share_number <- NULL
  total_uniq_ratio_to_them <- NULL
  total_uniq_ratio_to_us <- NULL
  total_share_ratio <- NULL
  total_sample <- NULL
  total_denomitor <- NULL
  for (samp in intercet_sample){
    
    our_each_mut <- our_MT[sample_name == samp, .(chrom_loc)]
    publish_each_mut <- publisMT[Case == samp,.(chrom_loc)]
    
    
    our_each <- nrow(our_MT[sample_name==samp,])
    sanger_each <- nrow(publisMT[Case==samp,])
    intercet_data <- nrow(intersect(our_each_mut,publish_each_mut))
    denominator <- our_each+sanger_each-intercet_data
    
    # count
    number_overlap <- nrow(intersect(our_each_mut,publish_each_mut))
    uniq_number_to_us <-nrow(setdiff(our_each_mut,publish_each_mut)) 
    uniq_number_to_them <- nrow(setdiff(publish_each_mut,our_each_mut)) 
    
    # ratio
    uniq_ratio_to_them <- nrow(setdiff(publish_each_mut,our_each_mut))/(denominator)
    uniq_ratio_to_us <- nrow(setdiff(our_each_mut,publish_each_mut))/(denominator)
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
    
  }
  
  data <- data.frame(share_ratio = as.numeric(total_share_ratio), 
                     sample = total_sample,
                     uniq_ratio_to_uga = as.numeric(total_uniq_ratio_to_us),
                     uniq_ratio_to_publication = as.numeric(total_uniq_ratio_to_them),
                     uniq_num_to_publication = as.numeric(total_uniq_num_to_them),
                     uniq_num_to_uga = as.numeric(total_uniq_num_to_us),
                     share_number = as.numeric(total_share_number),
                     total_denomitor= as.numeric(total_denomitor))
  data <- setDT(data)
  
  return (data)
}





## Anayzed used Mutect2 before mutation Distribution
base <- "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Compare_publication\\MT_mutateion_compare_with_korean"

publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep="\\"))
our_total <- before_total_summary

#publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
our_MT <- our_total
publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
our_MT <- clean_table(our_MT)
our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
our_MT$sample_name <- sapply(our_MT$sample_name, convert_MT_sample)
samples <- unique(publishMT$Case)


our_mut_number <- c()
pub_mut_number <- c()
sample_name <- c()

for (samp in samples){
  our_num <- nrow(our_total[sample_name==samp,]) 
  pb_num <- nrow(publishMT[Case==samp,])
  our_mut_number <- c(our_mut_number, our_num)
  pub_mut_number <- c(pub_mut_number,pb_num)
  sample_name <- c(sample_name,samp)
  
}

data <- data.frame(our_mut_count= our_mut_number,
                   pub_mut_number= pub_mut_number, 
                   sample_name = sample_name )



pdf(paste(base,"Before_filtering_Compare_MT_mut_mutect2.pdf",sep="\\")
    , height=4.84, width=4.84);

ggplot(data= data, aes(x=our_mut_count, y= pub_mut_number))+
  geom_point(shape=1,size= 3.5)+
  # geom_point(data=highlight_df, 
  #            aes(x=our_mut_count,y=sanger_mut_count), 
  #            color='red',
  #            size=3)+
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


## Anayzed used Mutect2 after mutation distribution

publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep="\\"))
our_total <- after_total_summary

#publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
our_MT <- our_total
publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
our_MT <- clean_table(our_MT)
our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
our_MT$sample_name <- sapply(our_MT$sample_name, convert_MT_sample)
samples <- unique(publishMT$Case)


our_mut_number <- c()
pub_mut_number <- c()
sample_name <- c()

for (samp in samples){
  our_num <- nrow(our_total[sample_name==samp,]) 
  pb_num <- nrow(publishMT[Case==samp,])
  our_mut_number <- c(our_mut_number, our_num)
  pub_mut_number <- c(pub_mut_number,pb_num)
  sample_name <- c(sample_name,samp)
  
}

data <- data.frame(our_mut_count= our_mut_number,
                   pub_mut_number= pub_mut_number, 
                   sample_name = sample_name )



pdf(paste(base,"After_filtering_Compare_MT_mut_mutect2.pdf",sep="\\")
    , height=4.84, width=4.84);

ggplot(data= data, aes(x=our_mut_count, y= pub_mut_number))+
  geom_point(shape=1,size= 3.5)+
  # geom_point(data=highlight_df, 
  #            aes(x=our_mut_count,y=sanger_mut_count), 
  #            color='red',
  #            size=3)+
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



## Analyzed Mut1 distribution

## after
publishMT <- fread(paste(base,"publish_MT_cds.txt",sep="\\"))
our_total <- fread(paste(base,"After_Total_VAF_depth_summary.txt",sep = "\\"));

publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
our_MT <- our_total[tumor_type=="MT",]


publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
our_MT <- clean_table(our_MT)
our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
our_MT$Case <- sapply(our_MT$sample_name, convert_MT_sample)
samples <- unique(publishMT$Case)


our_mut_number <- c()
pub_mut_number <- c()
sample_name <- c()

for (samp in samples){
  our_num <- nrow(our_total[sample_name==samp,]) 
  pb_num <- nrow(publishMT[Case==samp,])
  our_mut_number <- c(our_mut_number, our_num)
  pub_mut_number <- c(pub_mut_number,pb_num)
  sample_name <- c(sample_name,samp)
  
}

data <- data.frame(our_mut_count= our_mut_number,
                   pub_mut_number= pub_mut_number, 
                   sample_name = sample_name )



pdf(paste(base,"After_Compare_MT_mut.pdf",sep="\\")
    , height=4.84, width=4.84);

ggplot(data= data, aes(x=our_mut_count, y= pub_mut_number))+
  geom_point(shape=1,size= 3.5)+
  # geom_point(data=highlight_df, 
  #            aes(x=our_mut_count,y=sanger_mut_count), 
  #            color='red',
  #            size=3)+
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
## before mutect1
publishMT <- fread(paste(base,"publish_MT_cds.txt",sep="\\"))
our_total <- fread(paste(base,"Before_Total_VAF_sum.txt",sep = "\\"));

publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")
our_MT <- our_total[tumor_type=="MT",]


publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
our_MT <- clean_table(our_MT)
our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
our_MT$Case <- sapply(our_MT$sample_name, convert_MT_sample)
samples <- unique(publishMT$Case)


our_mut_number <- c()
pub_mut_number <- c()
sample_name <- c()

for (samp in samples){
  our_num <- nrow(our_total[sample_name==samp,]) 
  pb_num <- nrow(publishMT[Case==samp,])
  our_mut_number <- c(our_mut_number, our_num)
  pub_mut_number <- c(pub_mut_number,pb_num)
  sample_name <- c(sample_name,samp)
  
}

data <- data.frame(our_mut_count= our_mut_number,
                   pub_mut_number= pub_mut_number, 
                   sample_name = sample_name )



pdf(paste(base,"Before_Compare_MT_mut.pdf",sep="\\")
    , height=4.84, width=4.84);

ggplot(data= data, aes(x=our_mut_count, y= pub_mut_number))+
  geom_point(shape=1,size= 3.5)+
  # geom_point(data=highlight_df, 
  #            aes(x=our_mut_count,y=sanger_mut_count), 
  #            color='red',
  #            size=3)+
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


## Analyzed MT six bases before and after
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
source("C:\\Users\\abc73_000\\Documents\\GitHub\\Pan-Cancer\\six_base_function_util.R")

regular.text <- element_text(colour="black",size=20)
fontsize <- 20;
xangle <- 90
signature_colors <- c("cyan","black","red","gray","green","pink");
signature_levels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
base <- "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\VAF"

before_total <- fread(paste(base,"Before_Total_VAF_sum.txt",sep="\\"))
after_total <- fread(paste(base,"After_Total_VAF_sum.txt",sep="\\"))

before_MT <- before_total[tumor_type=="MT",]
after_MT <- after_total[tumor_type=="MT",]

total_sample <- sort(unique(before_MT$sample_name))

mut_before_count_sum <- NULL
mut_after_count_sum <- NULL

for (sample in total_sample){
  before <- conver_to_six_bases_basedon_mutation(before_MT,sample,"MT",rate=F)
  after <- conver_to_six_bases_basedon_mutation(after_MT,sample,"MT",rate = F)
  
  mut_before_count_sum <- rbind(mut_before_count_sum,before)
  mut_after_count_sum <- rbind(mut_after_count_sum,after)
}



mut_before_count_sum <- mut_before_count_sum[x!="CMT-33",]
mut_after_count_sum <- mut_after_count_sum[x!="CMT-33",]
## Mutect1 before

library(ggforce)
library(tidyr) 
library(dplyr)
library(ggplot2)
library(ggforce)

pdf(paste(base,"MT_samples_data_6bases_count.pdf",sep="\\")
    , height=16.94, width=12.94);

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
             nrow = 2)


dev.off()



# ### cat all of Mutect2 data together because the publication used Mutect2


cancer_type <- "MT"
base_dir = #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate/VAF/Mutect1"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Compare_publication/MT_mutateion_compare_with_korean/Mutect2/MT"

# # coverage_file <- read_excel(#"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate/New_WES_QC_dataset.xlsx",
# #   "C:\\Users\\abc73_000\\Desktop\\New_WES_QC_dataset.xlsx",
# #   sheet = "Mut")
#
# coverage_file <- setDT(coverage_file)

before_total_summary <- NULL
after_total_summary <- NULL

for (cancer in cancer_type){

  sample_names <- list.files(base_dir)

  for (sample in sample_names){

    M2before_file_name <- paste(sample,"vaf_before.txt",sep = '_')
    M2after_file_name <- paste(sample,"vaf_after.txt",sep = '_')

    Mutect2_before <-fread(paste(base_dir,sample,M2before_file_name ,sep=seperator))
    Mutect2_after <- fread(paste(base_dir,sample,M2after_file_name ,sep=seperator))
    # before
    if (nrow(Mutect2_before)!=0){
      before_each_file <-fread(paste(base_dir,sample,M2before_file_name ,sep=seperator))
      colnames(before_each_file) <- c("chrom","pos","vaf","ref","alt","tRef","tAlt","sample_name")
      before_each_file$sample_name <-sample
    }


    if (nrow(Mutect2_after)!=0){
      
      # after
      after_each_file <- fread(paste(base_dir,sample,M2after_file_name ,sep=seperator))
      colnames(after_each_file) <- c("chrom","pos","vaf","ref","alt","tRef","tAlt","sample_name")
      after_each_file$sample_name <-sample

    }
    before_total_summary <- rbind(before_total_summary,before_each_file)
    after_total_summary <- rbind(after_total_summary,after_each_file)
  }
}

before_total_summary$sample_name
## cat all file and write output

fwrite(before_total_summary,file = paste(base_dir,"MT_mutect2_before.txt",sep = seperator),
       sep="\t", col.names = T, row.names = F, quote = F)


fwrite(after_total_summary,file = paste(base_dir,"MT_mutect2_after.txt",sep = seperator),
       sep="\t", col.names = T, row.names = F, quote = F)


unique(before_total_summary$sample_name)
