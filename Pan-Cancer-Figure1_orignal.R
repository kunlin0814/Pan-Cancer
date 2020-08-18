library(tidyverse)
library(readxl)
library(wesanderson)
library(RColorBrewer)
#colors <- brewer.pal(3, "Set1"); # red, blue, green

## whole criteria 
#Sequence read pairs < 5M	
#Unique CDS mapping rate < 0.3	
#Mean_coverage < 30	
#Callable bases < 10Mb (samples/cases)

########### Photoshop use res = 450 #################




PAIR <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Figure1\\Original_Data_summary.xlsx",
                   sheet ="PAIRS")

callable <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Mutation_rate\\Mutation_rate.xlsx",
                       sheet ='All_samples_Retro_non_retro')


exclude <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Original_Data_summary.xlsx",
                      sheet ="Before_Matching_excluded")

total_file <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Methods_legends_tables\\Ploting_table.xlsx",
                         sheet ='NGS_Data_summary', skip =0)

#plot_result <- png("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\sequence_pairs.png",
#                   width=2800,height=1800,res=450)
regular.text <- element_text(colour="black",size=20);

tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\sequence_pairs.tiff", 
     width = 3500, height =1900, units = "px", res = 400)


pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Figure1.pdf"
    , height=4.8, width=6.2);

total_file %>% 
  #filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("MC","GLM","BCL","TCL","OM", "OSA","HSA","Unclassified")),
             y=as.numeric(Total_pairs)/1000000,fill=Status,color=Status)) +
  geom_point(size=1.8,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("Sequence Read Pairs in Millions")+
  scale_y_continuous(breaks = c(0,100,200))+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        legend.key=element_blank(),
        #legend.text =element_text(color = c("firebrick","black")),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_text(colour="black",size=22),
        #axis.title.y = element_text(colour="black",size=20,margin = margin(1,40,0,0)),
        text = element_text(colour="black",size=22),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour=c("black"),size=20,vjust=1.,hjust = .95, angle = 20),
        #axis.text.y = element_blank(),
        axis.text.y = element_text(colour="black",size=22),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("darkblue","red3"))+
  scale_fill_manual(values=c("darkblue","red3"))+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 19)+
  coord_cartesian(ylim=c(0,200))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=5, linetype="dashed", color = "yellow4", size = 1.5)

dev.off()

### Unique Mapping Rates ###
#plot_result <- png("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer_mapping_result\\Supplement_Figure1\\unique-mapping_pairs.png",
#                   width=2800,height=1800,res=350)

tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\unique-mapping_pairs.tiff", 
     width = 3500, height =1900, units = "px", res = 400)

total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("MC","GLM","BCL","TCL","OM", "OSA","HSA","Unclassified")),
             y=as.numeric(Uniquely_mapped_rate),fill=Status,color=Status)) +
  geom_point(size=1.8,position = position_jitterdodge(jitter.width = 0.3)) +
  ylab("Uniquely Mapping Rate")+
  scale_y_continuous(breaks = c(0,0.5,0.75,1.0))+
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
        axis.title.y = element_text(colour="black",size=24,margin = margin(1,5,.5,0)),
        text = element_text(colour="black",size=14),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour=c("black"),size=14,vjust=1,hjust = 0.9, angle = 45),
        axis.text.y = element_text(colour="black",size=24),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("darkblue","red3"))+
  scale_fill_manual(values=c("darkblue","red3"))+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 19)+
  coord_cartesian(ylim=c(0.3,1))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=0.6, linetype="dashed", color = "yellow4", size = 1.5)


dev.off()
###### gt 30 #####

tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\gt-30.tiff", 
     width = 3500, height =1900, units = "px", res = 400)

total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN)  %>% 
  #filter(!gt_30_fraction < 0.25 |gt_30_fraction ==NaN) %>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("MC","GLM","BCL","TCL","OM", "OSA","HSA","Unclassified")),
             y=as.numeric(gt_30_fraction),fill=Status,color=Status)) +
  geom_point(size=1.8,position = position_jitterdodge(jitter.width = 0.3)) +
  ylab("Fraction of mapping quality >30")+
  scale_y_continuous(breaks = c(0,0.5,0.75,1.0))+
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
        axis.title.y = element_text(colour="black",size=24,margin = margin(1,5,.5,0)),
        text = element_text(colour="black",size=14),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour=c("black"),size=14,vjust=1,hjust = 0.9, angle = 45),
        axis.text.y = element_text(colour="black",size=24),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("darkblue","red3"))+
  scale_fill_manual(values=c("darkblue","red3"))+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 19)+
  coord_cartesian(ylim=c(0.3,1))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))


dev.off()



#### Unique CDS ####
tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\unique-CDS-mapping.tiff", 
     width = 3500, height =1900, units = "px", res = 400)

total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(as.numeric(Uniquely_mapped_rate) >=0.6)%>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("MC","GLM","BCL","TCL","OM", "OSA","HSA","Unclassified")),
             y=as.numeric(uniq_CDS_region_paris_rates),fill=Status,color=Status)) +
  geom_point(size=1.8,position = position_jitterdodge(jitter.width = 0.3)) +
  
  ylab("Unique CDS Mapping Rates")+
  scale_y_continuous(breaks = c(0,0.4,0.8))+
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
        axis.title.y = element_text(colour="black",size=24,margin = margin(1,5,.5,0)),
        text = element_text(colour="black",size=14),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour=c("black"),size=14,vjust=1,hjust = 0.9, angle = 45),
        axis.text.y = element_text(colour="black",size=24),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("darkblue","red3"))+
  scale_fill_manual(values=c("darkblue","red3"))+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 19)+
  coord_cartesian(ylim=c(0,0.8))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=0.3, linetype="dashed", color = "yellow4", size = 1.5)

dev.off()

### Mean Coverage ###

tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\coverage-Mean.tiff"
     ,width = 3500, height =1900, units = "px", res = 400)

total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(as.numeric(Uniquely_mapped_rate) >=0.6) %>% 
  filter(as.numeric(uniq_CDS_region_paris_rates) >=0.3)%>%
  #filter(mean >30) %>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("MC","GLM","BCL","TCL","OM", "OSA","HSA","Unclassified")),
             y=as.numeric(mean),fill=Status,color=Status)) +
  geom_point(size=1.8,position = position_jitterdodge(jitter.width = 0.3)) +
  scale_y_continuous(breaks = c(0,100,200))+
  ylab("Mean Coverage")+
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
        axis.title.y = element_text(colour="black",size=24,margin = margin(1,5,.5,0)),
        text = element_text(colour="black",size=14),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour=c("black"),size=14,vjust=1,hjust = 0.9, angle = 45),
        axis.text.y = element_text(colour="black",size=24),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("darkblue","red3"))+
  scale_fill_manual(values=c("darkblue","red3"))+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 19)+
  coord_cartesian(ylim=c(0,200))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))+
  geom_hline(yintercept=30, linetype="dashed", color = "yellow4", size = 1.5)


dev.off()



#### RMSE ####

tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Coverage-RMSE.tiff", 
     width = 3500, height =1900, units = "px", res = 400)

total_file %>% 
  filter(!Total_pairs < 5000000 | Total_pairs==NaN) %>%
  filter(as.numeric(Uniquely_mapped_rate) >=0.6) %>% 
  filter(as.numeric(uniq_CDS_region_paris_rates) >=0.3)%>% 
  filter(mean>30) %>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("MC","GLM","BCL","TCL","OM", "OSA","HSA","Unclassified")),
             y=as.numeric(RMSE),fill=Status,color=Status)) +
  geom_point(size=1.8,position = position_jitterdodge(jitter.width = 0.3)) +
  scale_y_continuous(breaks = c(0,0.005,0.01))+
  ylab("RMSE of Sequence Coverage")+
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
        axis.title.y = element_text(colour="black",size=24,margin = margin(1,5,0,0)),
        text = element_text(colour="black",size=14),
        #axis.text.x = element_text(colour=c("black"),size=14,vjust=1,hjust = 0.9, angle = 45),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour="black",size=24),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.4, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("darkblue","red3"))+
  scale_fill_manual(values=c("darkblue","red3"))+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  scale_shape_manual(values = 19)+
  coord_cartesian(ylim=c(0,0.01))+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))


dev.off()


### change the color of the point ###
tiff(file = "G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\callable-bases.tiff", 
     width = 3500, height =1900, units = "px", res = 400)


total_file %>% 
  filter(Status=='Tumor') %>% 
  filter(!Case_ID %in% exclude$Cases) %>% 
  filter(!Callable_bases =='No-Pair')%>% 
  ggplot(aes(x=factor(Cancer_Type,levels = c("MC","GLM","BCL","TCL","OM", "OSA","HSA","Unclassified")),
             y=as.numeric(Callable_bases)/1000000,color='black')) +
  geom_point(size=1.8,position = position_jitterdodge(jitter.width = 0.3))+
  ylab("Callable bases in Millions")+
  scale_y_continuous(breaks = c(0,10,20,30))+
  theme(axis.line = element_line(colour = "black"),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key=element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.border = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(colour="black",size=24,margin = margin(1,1,0.5,0)),
        plot.margin = margin(1, 10, 4, 5),
        text = element_text(colour="black",size=20),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(colour=c("black"),size=14,vjust=1,hjust = 0.9, angle = 45),
        axis.text.y = element_text(colour="black",size=24),
        #axis.text.y = element_blank(),
        panel.background = element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.4, width = .7,colour = "black")+
  scale_shape_manual(values = 19)+
  scale_color_manual(values = 'black')+
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))+
  coord_cartesian(ylim=c(5,30))+
  geom_hline(yintercept=10, linetype="dashed", color = "yellow4", size = 1.5)
dev.off()
#scale_x_discrete(labels=panel2label)+
#scale_color_manual(values=wes_palette(n=5, name="Set3")) ## for scatter plot Wes Anderson color palettes
#scale_fill_manual(values=wes_palette(n=5, name="Royal2")) ## for box plot, bar plot Wes Anderson color palettes
#scale_shape_manual(values = 20)





