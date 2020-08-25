library(tidyverse)
library(data.table)
library(readxl)
library(gridExtra)

regular.text <- element_text(colour="black",size=20);
GeneCoverage <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/IntervalMeanCoverage/Interval_Log2_MDM2_CDKN2_DepthOfCoverage.xlsx",
                           sheet = "Log2Ratio_Interval")
GeneCoverage <- data.table(GeneCoverage)

geneName <- unique(GeneCoverage$GeneName)
sampleName <- unique(GeneCoverage$SampleName)
MDM2 <- c()
CDKN2A <- c()
CDKN2B <- c()
totalSummary <- NULL
for (sample in sampleName){
  MDM2log2 <- mean(GeneCoverage[(sampleName ==  sample) & (GeneName=="MDM2"),Log2Ratio])
  CDK2Alog2 <- mean(GeneCoverage[sampleName == sample & GeneName=="CDKN2A",Log2Ratio])
  CDK2Blog2 <- mean(GeneCoverage[sampleName == sample & GeneName=="CDKN2B",Log2Ratio])
  eachCol <- c(sample,MDM2log2, CDK2Alog2,CDK2Blog2)
  totalSummary <- rbind(totalSummary,eachCol)
}
totalSummary <- as.data.frame(totalSummary) 
# %>% 
#   write.table(file='C:\\Users\\abc73_000\\Desktop\\MDM2_mean.txt',
#               sep ='\t',row.names = F,quote = F, col.names = F)



GeneCoverage$SampleName

mean(GeneCoverage[SampleName == "DD0001" & GeneName=="MDM2",Log2Ratio])

GeneCoverage <- data.table(GeneCoverage)


#Median <- GeneCoverage[, .(MedianMDM2 = median(MDM2), MdeidanCDK2A = median(CDKN2A), MEdianCDKN2B = median(CDKN2B)),by = .(CancerType)]
#Mean <- GeneCoverage[, .(MeanMDM2 = mean(MDM2), MeanCDK2A = mean(CDKN2A), MeanCDKN2B = mean(CDKN2B)),by = .(CancerType)]


GeneCoverage[, .(.N),by = .(CancerType)]
Gene <- "MDM2"
GeneCoverage[with(GeneCoverage,Gene >1),]


GeneCoverage %>% 
  group_by(CancerType) %>% 
  summarise(Median = median(MDM2))

x <- factor(GeneCoverage$CancerType, levels=c( "OM", "OSA"));
fill <- factor(GeneCoverage$CancerType, levels=c( "OM", "OSA")); # how do you want to seperate the data, here use Noraml tumor to seperate
y <- as.numeric(GeneCoverage$MDM2)
data <- data.frame(x=x, y=y, fill=fill);
fill.colors <- c("#103B78", "#A0111A");
ylab <- "Log2 (Tumor/Normal)";



p <- get_jitter(data, x, y, fill, fill.colors=fill.colors, ylab=ylab, y_cutoffs = 10, 
                show.legend=F, xangle=30, compare_fills=FALSE, show.median=TRUE, dot_size=1.6);
p

pdf("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/IntervalMeanCoverage/Interval_OM_OSALog2Ration.pdf"
    , height=4.98, width=6.5);

GeneCoverage %>% 
  ggplot(aes(x=factor(CancerType,levels = c("OM", "OSA")),
             y=as.numeric(MDM2),color='black')) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("MDM2 Depth Coverage Log2 Ratio")+
  theme(axis.text=regular.text, 
        axis.title.y=element_text(colour="black",size=20, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.4, width = .7,colour = "black")+
  #scale_shape_manual(values = 19)+
  scale_color_manual(values = 'black')+
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))
GeneCoverage %>% 
  ggplot(aes(x=factor(CancerType,levels = c("OM", "OSA")),
             y=as.numeric(CDKN2A),color='black')) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("CDKN2A Depth Coverage Log2 Ratio")+
  theme(axis.text=regular.text, 
        axis.title.y=element_text(colour="black",size=20, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.4, width = .7,colour = "black")+
  #scale_shape_manual(values = 19)+
  scale_color_manual(values = 'black')+
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))
GeneCoverage %>% 
  ggplot(aes(x=factor(CancerType,levels = c("OM", "OSA")),
             y=as.numeric(CDKN2B),color='black')) +
  geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
  ylab("CDKN2B Depth Coverage Log2 Ratio")+
  theme(axis.text=regular.text, 
        axis.title.y=element_text(colour="black",size=20, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x =element_blank(),
        axis.text.x = element_text(angle=30, hjust=1), 
        panel.background=element_blank(), 
        axis.line=element_line(color="black"),
        legend.title=regular.text, 
        legend.position="none", 
        legend.text=regular.text, 
        legend.key=element_blank())+
  stat_summary(fun = median, fun.min = median, fun.max = median,position = "dodge",
               geom = "crossbar",size=0.4, width = .7,colour = "black")+
  #scale_shape_manual(values = 19)+
  scale_color_manual(values = 'black')+
  theme(plot.margin = unit(c(0.5,0.3,1,0.5), "cm"))
dev.off()





### Whole Region of Log2 Ratio ###

file <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/PanCancerCoverage/MC_Coverage/CMT-13_CoverageLog2Ratio.txt",
              header = F)
chrom_pos <- strsplit(file$V1,split = ":")

categorize <- function(x){
 data <- strsplit(x,split = ":")
 return (data[[1]][1])
}

chromCategory <- as.vector(sapply(file$V1,FUN =categorize ))

file$V3 <- chromCategory
chrom <- c()
for (i in 1:38){
  chrom <- c(chrom,paste("chr",i,sep = ""))
}
chrom <- c(chrom, "chrX")

regular.text <- element_text(colour="black",size=20);

pdf("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/PanCancerCoverage/MC_Coverage/V6MC_Log2Ratio.plot.pdf"
    , height=4.98, width=6.5);
testData <- file[1:10000]

p <- ggplot(data=file, aes(x=factor(file$V1, levels = file$V1), y=V2, fill = V3)) +
  geom_bar(stat="identity")+
  xlab("Chrom Position")+
  ylab("log2 Ratio (Tumor/Normal)")+
  ggtitle("CMT-2")+
  theme_minimal()+
  theme(
    legend.title= element_blank(), 
    legend.position="top", 
    legend.text=element_text(size = 10), 
    title = regular.text,
    axis.title.x = regular.text,
    axis.title.y = regular.text,
    axis.text.x = element_blank(),
    axis.text.y = regular.text,
    panel.background=element_blank(),
    legend.key.size = unit(0.25,"line"),
    axis.line=element_line(color="black"))
print(p)


dev.off()




plotCoverage <- function(file1,file2,NumberDataPoint){
  Normal <- fread(file1,header= F)
  Tumor <- fread(file2,header= F)
  
  
  Normal$V1 <- factor(Normal$V1)
  Tumor$V1 <- factor(Tumor$V1)
  
  new_Normal <- data.frame(Normal[1:NumberDataPoint])
  new_Tumor <- data.frame(Tumor[1:NumberDataPoint])
  
  
  
  
  p1 <- ggplot(data=new_Normal, aes(x=V1, y=V2)) +
    geom_bar(stat="identity", fill="steelblue")+
    xlab("Chrom Position")+
    ylab("Mean Coverage")+
    ggtitle(file1)+
    theme(
      axis.title.x =element_blank(),
      axis.text.x = element_blank())+
    coord_cartesian(ylim=c(0,500))+
    theme_minimal()
  
  
  p2 <- ggplot(data=new_Tumor, aes(x=V1, y=V2)) +
    geom_bar(stat="identity", fill="steelblue")+
    xlab("Chrom Position")+
    ylab("Mean Coverage")+
    ggtitle(file2)+
    theme(
      axis.title.x =element_blank(),
      axis.text.x = element_blank())+
    coord_cartesian(ylim=c(0,500))+
    theme_minimal()
  
  
  grid.arrange(p1, p2, ncol=1, nrow= 2)  
}