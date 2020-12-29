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
dot_size <- 1.4;
abs_text_size <- 16;
regular.text <- element_text(colour="black",size=abs_text_size);
xangle <- 45;
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
regular.text <- element_text(colour="black",size=20)

base <- "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate_VAF\\Mut_rate\\HSA_issue"

HSA_data <- read_excel(paste(base,"Kun-Lin_Broad_HSAs_Mutations.xlsx",sep = "\\"),sheet = "mut_rate_compare")
HSA_data <- setDT(HSA_data)
final_HSA <- HSA_data[!Case %in% excldue$Cases & Case !="HSA_9",]
samples <- final_HSA$Case


pdf(paste(base,"after_filtering_Compare_HSA_mut_mutect12.pdf",sep="\\")
    , height=4.84, width=4.84);

ggplot(data= final_HSA, aes(x=as.numeric(our_data_mutect1), y= as.numeric(Broad_data)))+
  geom_point(shape=1,size= 3.5)+
  # geom_point(data=highlight_df, 
  #            aes(x=our_mut_count,y=sanger_mut_count), 
  #            color='red',
  #            size=3)+
  ylim(0,125)+
  xlim(0,125)+
  xlab("Our Mutect1")+
  ylab("Broad Mutect2")+
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




ggplot(data= final_HSA, aes(x=as.numeric(our_mutect2_5steps), y= as.numeric(Broad_data)))+
  geom_point(shape=1,size= 3.5)+
  # geom_point(data=highlight_df, 
  #            aes(x=our_mut_count,y=sanger_mut_count), 
  #            color='red',
  #            size=3)+
  ylim(0,125)+
  xlim(0,125)+
  xlab("Our_mutect2")+
  ylab("Broad_mutect2")+
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






## Analyze MT samples

base <- "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Compare_publication\\MT_mutateion_compare_with_korean"
publishMT <- fread(paste(base,"PON_DbSNP_CDS_combineChrompublishMT.txt",sep="\\"))
our_total <- fread(paste(base,"After_Total_VAF_depth_summary.txt",sep = "\\"));
#publishMT$Chromosome <- paste("chr",publishMT$Chromosome,sep = "")


our_MT <- our_total[tumor_type=="MT",]


publishMT$chrom_loc <- paste(publishMT$Chromosome,publishMT$Position,sep = "_")
publisMT <- setDT(publishMT)
our_MT <- clean_table(our_MT)
our_MT$chrom_loc <- paste(our_MT$chrom,our_MT$pos,sep= "_")
our_MT$sample_name <- sapply(our_MT$sample_name, convert_MT_sample)
samples <- unique(publishMT$Case)
our_sample <- unique(our_MT$sample_name)

their_sample <- unique(publishMT$Case)
intercet_sample <- intersect(their_sample,our_sample)
#setdiff(their_sample,our_sample)

## Anayzed ratio
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
