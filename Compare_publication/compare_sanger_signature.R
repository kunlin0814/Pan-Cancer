library(tidyverse)
library(data.table)
library(gridExtra)
library(ggrepel)
library(readxl)
library(ggpubr)
library(grid)

source("C:\\Users\\abc73_000\\Documents\\GitHub\\VAF\\six_base_function_util.R")

excldue <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan-Cancer-Manuscript\\Figure1\\Original_Data_summary.xlsx",
                      sheet = "Before_Matching_excluded")
fontsize=20
dot_size <- 1.4;
abs_text_size <- 16;
regular.text <- element_text(colour="black",size=abs_text_size);
xangle <- 45;
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
regular.text <- element_text(colour="black",size=20)


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

#### Analyzed the ratio that overlap

base <- "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Compare_publication\\OM_mutation_compare_with_Sanger"
# sanger_signature <- read_excel(paste(base,"Sanger_mutation.xlsx",sep ="\\"),
#                                 sheet ="Canine", skip = 29)

sanger_signature <- fread(paste(base,"DbSNP_sanger_mut_file_before_DBSNP.txt",sep ="\\"))
clean_sanger <- sanger_signature
#clean_sanger <- sanger_signature[,c("#Chr","Position","Ref","Alt","Sample")]
#clean_sanger$Chromosome <- paste("chr",clean_sanger$`#Chr`,sep ="")
clean_sanger$chrom_loc <- paste(clean_sanger$Chromosome,clean_sanger$Position,sep = "_")
clean_sanger <- setDT(clean_sanger)

samples <- sort(unique(sanger_signature$Sample))
our_data <- fread(paste(base,"our_totalSangerOM_mutation.txt",sep="\\"))
colnames(our_data) <- c("Chromosome","Position","Ref","Alt","Sample")

our_data_sample_convert <- sapply(our_data$Sample,convert_sample)
our_data$Sample <- our_data_sample_convert
our_data$chrom_loc <- paste(our_data$Chromosome,our_data$Position,sep ="_")

total_uniq_num_to_them <- NULL
total_uniq_num_to_us <- NULL
total_share_number <- NULL
total_uniq_ratio_to_them <- NULL
total_uniq_ratio_to_us <- NULL
total_share_ratio <- NULL
total_sample <- NULL
total_denomitor <- NULL
for (samp in samples){

our_each_mut <- our_data[Sample == samp, .(chrom_loc)]
sanger_each_mut <- clean_sanger[Sample == samp,.(chrom_loc)]


our_each <- nrow(our_data[Sample==samp,])
sanger_each <- nrow(clean_sanger[Sample==samp,])
intercet_data <- nrow(intersect(our_each_mut,sanger_each_mut))
denominator <- our_each+sanger_each-intercet_data

# count
number_overlap <- nrow(intersect(our_each_mut,sanger_each_mut))
uniq_number_to_us <-nrow(setdiff(our_each_mut,sanger_each_mut)) 
uniq_number_to_them <- nrow(setdiff(sanger_each_mut,our_each_mut)) 

# ratio
uniq_ratio_to_them <- nrow(setdiff(sanger_each_mut,our_each_mut))/(denominator)
uniq_ratio_to_us <- nrow(setdiff(our_each_mut,sanger_each_mut))/(denominator)
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
                   uniq_ratio_to_sanger = as.numeric(total_uniq_ratio_to_them),
                   uniq_num_to_sanger = as.numeric(total_uniq_num_to_them),
                   uniq_num_to_uga = as.numeric(total_uniq_num_to_us),
                   share_number = as.numeric(total_share_number),
                   total_denomitor= as.numeric(total_denomitor))
data <- setDT(data)

pdf(paste(base,"bar_OM_Mutation_overlap_ratio.pdf",sep="\\")
    , height=12.94, width=12.94);


count_data <- melt(data, id.vars = c("sample"),
                   measure.vars= c("uniq_num_to_uga","uniq_num_to_sanger","share_number"),
                   variable.name = "fill")
count_data <- count_data[order(sample)]

x <- count_data$sample
y <- count_data$value
classify <- c("uniq_num_to_uga","uniq_num_to_sanger","share_number");
fill <- count_data$fill
fill <- factor(fill, levels=classify);
samples <- unique(x);
sample_order <- order(sapply(samples, function(s) {sum(y[which(x == s)])}), decreasing=T);
x <- factor(x, levels=samples[sample_order]);
plot_data <- data.frame(x=x, y=y, fill=fill);
fill_colors <- c("cyan","black","red");

print(p)


ratio_data <- melt(data, id.vars = c("sample"),
                   measure.vars= c("uniq_ratio_to_uga","uniq_ratio_to_sanger","share_ratio"),
                   variable.name = "fill")
ratio_data <- ratio_data[order(sample)]

x <- ratio_data$sample
y <- ratio_data$value
classify <- c("uniq_ratio_to_uga","uniq_ratio_to_sanger","share_ratio");
fill <- ratio_data$fill
fill <- factor(fill, levels=classify);
samples <- unique(x);
sample_order <-sample_order;
x <- factor(x, levels=samples[sample_order]);
plot_data <- data.frame(x=x, y=y, fill=fill);
fill_colors <- c("cyan","black","red");

p <- ggplot(plot_data, aes(x=x, y=y, fill=fill)) + 
  geom_bar(stat="identity",position='stack', width=0.6)+
  ggtitle("OM mutation ratio overlapped with Sanger")+
  scale_fill_manual(values=fill_colors)+
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title.x = element_blank(),
    #element_text(face="plain",colour="black",size=fontsize),
    axis.title.y = element_blank(),
    #element_text(face="plain",colour="black",size=fontsize),
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size=fontsize,face="plain",colour="black"),
    legend.position="bottom",
    legend.title= element_blank(), legend.text = element_text(size=fontsize,face="plain",colour="black"))

print(p)

dev.off()
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
