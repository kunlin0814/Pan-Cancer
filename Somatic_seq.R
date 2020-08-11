library(ggplot2)
library(data.table)
library(readxl)

file <- fread("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\Mutect2_somatic_seq_\\Somatic_seq_total.txt")

exclude <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Original_Data_summary.xlsx",
                      sheet ="Before_Matching_excluded")



file[, .(Callable_bases)]




file <- file[!file_name %in% exclude$Cases ]
#Mammary <- file[Cancer_Type=='Mammary_Cancer']
regular.text <- element_text(colour="black",size=20);
pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Somatic_seq.pdf"
    , height=4.8, width=6.2);



ggplot(data= file,aes(x= Cancer_Type, y = (file$`Mut_PASS-1`)*1000000/(file$Callable_bases), color=Cancer_Type))+
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.1))+
  ylab("Mutation Rates")+
  ggtitle("Three Overlap")+
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        #axis.title.x =element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")


ggplot(data= file,aes(x= Cancer_Type, y = ((file$`Mut_PASS-2-1`)+(file$`Mut_PASS-2-2`)+(file$`Mut_PASS-3`))*100000/(file$Callable_bases),
                                                     color=Cancer_Type))+
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.1))+
  ylab("Mutation Rates")+
  ggtitle("Any Two")+
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        #axis.title.x =element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")


ggplot(data= file,aes(x= Cancer_Type, y = ((file$`Mut_PASS-2`)*100000)/(file$Callable_bases), color=Cancer_Type))+
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.1))+
  ylab("Mutation Rates")+
  ggtitle("Mutect2 Only")+
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        #axis.title.x =element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")

ggplot(data= file,aes(x= Cancer_Type, y = (file$`Mut_PASS-2`+file$`Mut_PASS-3`)*1000000/(file$Callable_bases),
                      color=Cancer_Type))+
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.1))+
  ylab("Mutation Rates")+
  ggtitle("Mutect2 Only and VC Overlap")+
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        #axis.title.x =element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  coord_cartesian(ylim=c(0,200))
dev.off()


