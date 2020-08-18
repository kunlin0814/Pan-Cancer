library(ggplot2)
library(data.table)
library(readxl)
library(gridExtra)

Create_VAF_plot <- function(before5steps,after5steps, title){
yscale <- max(before5steps$V3,after5steps$V3)
# if (max(before5steps$V3)> max(after5steps$V3)){
#   yscale <- max(before5steps$V3)
# }
# else{
#   yscale <- max(after5steps$V3)
# }
regular.text <- element_text(colour="black",size=14);
  
p1 <- ggplot(data= before5steps,aes(x=0, y = as.numeric(before5steps$V3),color='black'))+
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.05))+
  ylab("VAF")+
  ylim(0,yscale)+
  ggtitle("Before 5 steps Mutect")+
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
   stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
                geom = "crossbar",size=0.5, width = .7,colour = "black")

p2 <- ggplot(data= after5steps,aes(x=0, y = as.numeric(after5steps$V3),color='black'))+
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.05))+
  ylab("VAF")+
  ylim(0,yscale)+
  ggtitle("After 5 steps Mutect")+
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")

grid.arrange(p1, p2, nrow = 1, top = title)
}

### Mutect 2 ####
file <- list.files (path = "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate/VAF/Mutect2/Unclassified/")

# pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Unclassified_VAF_Mutect2.pdf"
#     , height=4.8, width=6.2);

for (i in file){
  path <- paste("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate/VAF/Mutect2/Unclassified/",i,sep="")
  before_file <- paste(i,"_VAF_Before.txt",sep="")
  after_file <- paste(i,"_VAF_After.txt",sep="")
  if (file.size(paste(path,before_file,sep="/"))> 0 & file.size(paste(path,after_file,sep="/")) > 0){
  before5steps <- fread(paste(path,before_file,sep="/"))
  after5steps <- fread(paste(path,after_file,sep="/"))
  Create_VAF_plot(before5steps,after5steps,i)
  }
  
}


### Mutect 1 ####
file <- list.files (path = "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate/VAF/Mutect1/TCELL/")

pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\TCELL_VAF_Mutect1.pdf"
    , height=4.8, width=6.2);

for (i in file){
  path <- paste("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate/VAF/Mutect1/TCELL/",i,sep="")
  before_file <- paste(i,"_VAF_Before.txt",sep="")
  after_file <- paste(i,"_VAF_After.txt",sep="")
  
  #print(before_file)
 
  #print(paste(path,before_file,sep="/"))
  #print(paste(path,after_file,sep="/"))
  if (file.size(paste(path,before_file,sep="/"))> 0 & file.size(paste(path,after_file,sep="/")) > 0){
    before5steps <- fread(paste(path,before_file,sep="/"))
    after5steps <- fread(paste(path,after_file,sep="/"))
    Create_VAF_plot(before5steps,after5steps,i)
  }
  
}


dev.off()


Mutect1_TotalVAF <- fread("C:\\Users\\abc73_000\\Desktop\\VAF\\Mutect1\\Total_Mutect1VAF.txt")
Mutect2_TotalVAF <- fread("C:\\Users\\abc73_000\\Desktop\\VAF\\Mutect2\\Total_Mutect2VAF.txt")

library(tidyverse)

pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Mutect12_VAF.pdf"
    , height=4.8, width=6.2);

p1 <- Mutect2_TotalVAF %>% 
  filter(V4=='Before') %>% 
  filter(V5=="BCELL") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect2 BCL Before VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)


p2 <- Mutect2_TotalVAF %>% 
  filter(V4=='After') %>% 
  filter(V5=="BCELL") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect2 BCL After VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)

p3 <-Mutect1_TotalVAF %>% 
  filter(V4=='Before') %>% 
  filter(V5=="BCELL") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect1 BCL Before VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)


p4 <-Mutect1_TotalVAF %>% 
  filter(V4=='After') %>% 
  filter(V5=="BCELL") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect1 BCL After VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)
grid.arrange(p1, p2, p3,p4,nrow = 2)



p1 <-Mutect2_TotalVAF %>% 
  filter(V4=='Before') %>% 
  filter(V5=="TCELL") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect2 TCL Before VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)
p2 <-Mutect2_TotalVAF %>% 
  filter(V4=='After') %>% 
  filter(V5=="TCELL") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect2 TCL After VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)

p3 <-Mutect1_TotalVAF %>% 
  filter(V4=='Before') %>% 
  filter(V5=="TCELL") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect1 TCL Before VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)
p4 <-Mutect1_TotalVAF %>% 
  filter(V4=='After') %>% 
  filter(V5=="TCELL") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect1 TCL After VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)


grid.arrange(p1, p2, p3,p4,nrow = 2)


p1 <- Mutect2_TotalVAF %>% 
  filter(V4=='Before') %>% 
  filter(V5=="MC") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect2 MC Before VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)

p2 <- Mutect2_TotalVAF %>% 
  filter(V4=='After') %>% 
  filter(V5=="MC") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect2 MC After VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)


p3 <-Mutect1_TotalVAF %>% 
  filter(V4=='Before') %>% 
  filter(V5=="MC") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect1 MC Before VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)

p4 <-Mutect1_TotalVAF %>% 
  filter(V4=='After') %>% 
  filter(V5=="MC") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect1 MC After VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)

grid.arrange(p1, p2, p3,p4, nrow = 2)

p1 <- Mutect2_TotalVAF %>% 
  filter(V4=='Before') %>% 
  filter(V5=="OM") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect2 OM Before VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)
p2 <- Mutect2_TotalVAF %>% 
  filter(V4=='After') %>% 
  filter(V5=="OM") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect2 OM After VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)
p3 <-Mutect1_TotalVAF %>% 
  filter(V4=='Before') %>% 
  filter(V5=="OM") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect1 OM Before VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)
p4 <-Mutect1_TotalVAF %>% 
  filter(V4=='After') %>% 
  filter(V5=="OM") %>% 
  ggplot(aes(x=V3)) + 
  geom_histogram(bins=5000)+
  xlab("VAF")+
  ggtitle("Mutect1 OM After VAF")+
  geom_vline(aes(xintercept = median(V3)),col='blue',size=1)

grid.arrange(p1, p2, p3,p4, nrow = 2)
dev.off()









Mutect1_TotalVAF[V4 == 'Before' ,.(.N),by=.(V5)]
Mutect1_TotalVAF[V4 == 'After' ,.(.N),by=.(V5)]

Mutect2_TotalVAF[V4 == 'Before' ,.(.N),by=.(V5)]
Mutect2_TotalVAF[V4 == 'After' ,.(.N),by=.(V5)]





pdf("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Mutect12_Before_AfterVAF.pdf"
    , height=4.8, width=6.2);
regular.text <- element_text(colour="black",size=14);
ggplot(data= Mutect1_TotalVAF, 
       aes(x=factor(Mutect1_TotalVAF$V5,levels = c("MC","BCELL","TCELL","OM")),
       y=as.numeric(Mutect1_TotalVAF$V3),fill=V4,color=V4)) +
  
  geom_point(size=3,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("VAF")+
  ggtitle("Mutect1 VAF")+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="top", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("darkblue","red3"))+
  scale_fill_manual(values=c("darkblue","red3"))+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 20)+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))
  #coord_cartesian(ylim=c(0,200))+


ggplot(data= Mutect2_TotalVAF, 
       aes(x=factor(Mutect2_TotalVAF$V5,levels = c("MC","BCELL","TCELL","OM")),
           y=as.numeric(Mutect2_TotalVAF$V3),fill=V4,color=V4)) +
  
  geom_point(size=3,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("VAF")+
  ggtitle("Mutect2 VAF")+
  #labs(subtitle = "p<0.01")+
  #stat_compare_means(aes(group=gene),label="p.signif",symnum.args = symnumargs,label.y = c(2.1,2.1,2.1,2.1,3,3)) +
  #scale_y_log10(limits=c(0.001,1000),breaks=c(0.001,0.01,0.1,1,10,100,1000),labels=c(0,0.01,0.1,1,10,100,1000)) +
  theme(axis.text=regular.text, 
        axis.title.y=regular.text,
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=0, hjust=0.5), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="top", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,position = "dodge",
               geom = "crossbar",size=0.5, width = .7,colour = "black")+
  #scale_x_discrete(labels=panel2label)+
  scale_color_manual(values = c("darkblue","red3"))+
  scale_fill_manual(values=c("darkblue","red3"))+
  #scale_fill_manual(values=c("firebrick","darkolivegreen"))+
  #scale_shape_manual(values = 20)+
  theme(plot.margin = unit(c(1,0.3,1.5,0.5), "cm"))
#coord_cartesian(ylim=c(0,200))+

dev.off()
Mutect1VAF_summary <- Mutect1_TotalVAF[, .(.N), by=.(V4,V5)]
Mutect2VAF_summary <- Mutect2_TotalVAF[, .(.N), by=.(V4,V5)]

total_reason_filter <- fread("C:\\Users\\abc73_000\\Desktop\\Total_Reasons_filteringout.txt")



