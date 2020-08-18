library(tidyverse)
library(data.table)
library(readxl)
library(gridExtra)

file <- fread("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/PanCancerCoverage/MC_Coverage/CMT-2/testlog2.txt",
              header = F)

pdf("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/PanCancerCoverage/MC_Coverage/MC_Log2Ratio.plot.pdf"
    , height=4.98, width=6.5);

ggplot(data=file[1:100000, ], aes(x=V1[1:100000], y=V2[1:100000])) +
  geom_bar(stat="identity", fill="steelblue")+
  xlab("Chrom Position")+
  ylab("log2 Ratio")+
  ggtitle("test")+
  theme(
    axis.title.x =element_blank(),
    axis.text.x = element_blank())+
  #coord_cartesian(ylim=c(0,500))+
  theme_minimal()

dev.off()


file[V2<100,]



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