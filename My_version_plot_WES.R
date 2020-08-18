library(tidyverse)
library(readxl)
library(wesanderson)
library(RColorBrewer)

total_file <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan-Cancer-Manuscript/Methods_legends_tables/V6Combined_Summary_Public_Data.xlsx",
                         sheet ='NGS_Data_summary', skip = 3)


exclude <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Original_Data_summary.xlsx",
                      sheet ="Before_Matching_excluded")

averge_randomness <- read_excel("C:\\Users\\abc73_000\\Desktop\\Validation Data_Set\\Validation_Set_Mutation.xlsx",
                                sheet ='Randomness')

gt30 <-read_excel("C:\\Users\\abc73_000\\Desktop\\Validation Data_Set\\Validation_Set_Mutation.xlsx",
                  sheet ='gt30')

callable <-  read_excel("C:\\Users\\abc73_000\\Desktop\\Validation Data_Set\\Validation_Set_Mutation.xlsx",
                        sheet ='callable')

sample_pair <- read_excel("C:\\Users\\abc73_000\\Desktop\\New_Data_OM_HSA_OSA\\NEW_DATA_QC.xlsx",
                          sheet ='Sample_pair')

HSA <-  read_excel("C:\\Users\\abc73_000\\Desktop\\New_Data_OM_HSA_OSA\\NEW_DATA_QC.xlsx",
                   sheet ='HSA')



#match("SRR6277192", sample_pair[1,], nomatch = 0)

check_sample_name <- function(x){
  total_row <- nrow(sample_pair)
  for (i in 1:total_row){
    if (x %in% sample_pair[i,]){
      return (sample_pair$Cases[i])
    }
  }
  return ("NaN")
}

HSA_sample <- sapply(as.vector(HSA$File_name),check_sample_name )
HSA$Case <- HSA_sample

HSA %>% 
  filter(Uniq_mapped_rate < 0.6) %>% 
  write.table(file='C:\\Users\\abc73_000\\Desktop\\HSAlt0.6unmapped',
              sep ='\t',row.names = F,quote = F, col.names = T)

HSA %>% 
  filter(Uniq_mapped_rate >= 0.6) %>%
  filter(mean <30) %>% 
  write.table(file='C:\\Users\\abc73_000\\Desktop\\meanlt30',
              sep ='\t',row.names = F,quote = F, col.names = T)


total_file <- unique(total_file) 

total_file %>% 
  filter(Uniq_mapped_rate <0.6)



exclude <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Original_Data_summary.xlsx",
                      sheet ="Before_Matching_excluded")






regular.text <- element_text(colour="black",size=24);
##### My version #######

tiff(file = "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Figure1/Supplementary/V6/Total_Sequence_pair.tiff",
     width = 3500, height =3000, units = "px", res = 400)

total_file %>% 
  #filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","LYM","OM", "OSA","HSA","UCL")),
             y=as.numeric(Total_pairs)/1000000,fill=Status,color=Status)) +
  
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Sequence Read Pairs in Millions")+
  #scale_y_continuous(breaks = c(0,100,200))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=element_blank(), 
        legend.position="none", 
        legend.text=element_blank(), 
        legend.key=element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("darkblue","red3"))+
  scale_fill_manual(values=c("darkblue","red3"))+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 20)+
  coord_cartesian(ylim=c(0,200))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=5, linetype="longdash", color = "yellow4", size = 0.7)

dev.off()



tiff(file = "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Figure1/Supplementary/V6/unique-mapping_pairs.tiff",
     width = 3500, height =3000, units = "px", res = 400)

total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","LYM","OM", "OSA","HSA","UCL")),
             y=as.numeric(Uniquely_mapped_rate),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Uniquely Mapping Rate")+
  scale_y_continuous(breaks = c(0,0.5,0.75,1.0))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=element_blank(), 
        legend.position="none", 
        legend.text=element_blank(), 
        legend.key=element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("darkblue","red3"))+
  scale_fill_manual(values=c("darkblue","red3"))+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 19)+
  coord_cartesian(ylim=c(0.3,1))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=0.6, linetype="longdash", color = "yellow4", size = 0.7)


dev.off()
###### gt 30 #####

tiff(file = "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Figure1/Supplementary/V6/gt-30.tiff",
     width = 3500, height =3000, units = "px", res = 400)

total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN)  %>% 
  #filter(!gt_30_fraction < 0.25 |gt_30_fraction ==NaN) %>% 
  ggplot(aes(x=factor(Tumor_Type,levels = c("MT","GLM","LYM","OM", "OSA","HSA","UCL")),
             y=as.numeric(gt_30_fraction),fill=Status,color=Status)) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Fraction of mapping quality >30")+
  scale_y_continuous(breaks = c(0,0.5,0.75,1.0))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=element_blank(), 
        legend.position="none", 
        legend.text=element_blank(), 
        legend.key=element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("darkblue","red3"))+
  scale_fill_manual(values=c("darkblue","red3"))+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 19)+
  coord_cartesian(ylim=c(0.3,1))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))


dev.off()

