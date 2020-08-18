library(tidyverse)
library(readxl)
library(wesanderson)
library(RColorBrewer)

Mut_data <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan_cancer_mapping_result/Mutation_rate/Mutation_rate.xlsx",
                       sheet ='All_samples_Retro+indel' )

total_file <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\Data_summary.xlsx",
                         sheet ='Total')

PAIR <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\Data_summary.xlsx",
                   sheet ='PAIRS')

callable <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\Data_summary.xlsx",
                       sheet ='non-retro-mut')

Lymphoma_subtype <- read.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Meta\\Pancancer_metadata_05_11_2020.txt",
                               sep ='\t',header =T, stringsAsFactors =F)

Final_table <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\Summary_of_public_data.xlsx",
                          sheet ='NGS_Data_summary')


######

gt_0.8 <- total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(Cancer_Type =='Glioma') %>% 
  filter(as.numeric(Uniquely_mapped_rate) >=0.8) %>% 
  select(ID)
#write.table("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Base_Quality\\Glioma_data\\gt0.8_glioma.txt",
#            sep ='\t',row.names = F,quote = F)

lt_0.8 <- total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(Cancer_Type =='Glioma') %>% 
  filter(as.numeric(Uniquely_mapped_rate) <0.8) %>%
  select(ID)


lt_0.8Cases <- PAIR %>% 
  filter(Cancer_type=='Glioma') %>% 
  filter(Normal %in% lt_0.8$ID | Tumor %in% lt_0.8$ID) %>% 
  select(Cases)

gt_0.8Cases <- PAIR %>% 
  filter(Cancer_type=='Glioma') %>% 
  filter(!(Normal %in% lt_0.8$ID | Tumor %in% lt_0.8$ID)) %>% 
  select(Cases)

gt0.8Mut_rate <- Mut_data %>% 
  filter(Cancer_Type=="Glioma") %>% 
  filter(Status == 'Non-Retro') %>%
  filter(file_name %in% gt_0.8Cases$Cases) %>% 
  mutate(Mapping = 'gt0.8')

lt0.8Mut_rate <- Mut_data %>% 
  filter(Cancer_Type=="Glioma") %>% 
  filter(Status == 'Non-Retro') %>%
  filter(file_name %in% lt_0.8Cases$Cases) %>% 
  mutate(Mapping ='lt0.8')

total_glioma_class_mut <- rbind(gt0.8Mut_rate,lt0.8Mut_rate)

png(file = "Mutation_Rate_Between_unique_mapping.tiff", width = 1700, height =1000, units = "px", res = 300)
total_glioma_class_mut %>% 
  #filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Mapping,levels = c("gt0.8","lt0.8")),
             y=as.numeric(Mutation_Rate))) +
  geom_point(size=0.8,aes(color=Mapping),position = position_jitterdodge(jitter.width = 0.1))+
  ylab("TMB")+
  ylim(0, max(total_glioma_class_mut$Mutation_Rate))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key=element_blank(),
        #legend.text =element_text(color = c("firebrick","black")),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=26,margin = margin(1,5,0.5,0)),
        plot.margin = margin(1, 10, 4, 5),
        text = element_text(colour="black",size=20),
        axis.text.x = element_text(colour=c("black"),size=14,vjust=1,hjust = 0.9, angle = 45),
        axis.text.y = element_text(colour="black",size=20),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  scale_color_brewer(palette="Dark2")+### for scatter plot
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))+
  coord_cartesian(ylim=c(0,2.5))

dev.off()


png(file = "Callable_Between_unique_mapping.tiff", width = 1700, height =1000, units = "px", res = 300)
total_glioma_class_mut %>% 
  #filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Mapping,levels = c("gt0.8","lt0.8")),
             y=as.numeric(Callable_bases))) +
  geom_point(size=0.8,aes(color=Mapping),position = position_jitterdodge(jitter.width = 0.1))+
  ylab("TMB")+
  ylim(0, max(total_glioma_class_mut$Callable_bases))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key=element_blank(),
        #legend.text =element_text(color = c("firebrick","black")),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=26,margin = margin(1,5,0.5,0)),
        plot.margin = margin(1, 10, 4, 5),
        text = element_text(colour="black",size=20),
        axis.text.x = element_text(colour=c("black"),size=14,vjust=1,hjust = 0.9, angle = 45),
        axis.text.y = element_text(colour="black",size=20),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.1, width = .8,colour = "black")+
  scale_color_brewer(palette="Dark2")+### for scatter plot
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))
  coord_cartesian(ylim=c(0,2.5))

dev.off()


wilcox.test(gt0.8Mut_rate$Mutation_Rate, lt0.8Mut_rate$Mutation_Rate)
#Pvalue = 0.15

median(gt0.8Mut_rate$Mutation_Rate)
median(lt0.8Mut_rate$Mutation_Rate)
