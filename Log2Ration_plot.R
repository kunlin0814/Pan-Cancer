library(tidyverse)
library(data.table)
library(readxl)
library(gridExtra)

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