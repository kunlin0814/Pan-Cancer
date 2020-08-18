library(tidyverse)
library(data.table)
library(readxl)
library(gridExtra)

#"C:\\Users\\abc73_000\\Desktop\\Bioproject_check"
#"G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Grant_table\\Cancer_Sample\\MC\\WXS\\"
pdf("C:/Users/abc73_000/Desktop/GLM_Coverage/GLM_Coverage.plot.pdf"
    , height=4.98, width=6.5);
location <- "C:/Users/abc73_000/Desktop/GLM_Coverage/"
Samples <- list.dirs("C:/Users/abc73_000/Desktop/GLM_Coverage/")
length(Samples)
for (i in 2:length(Samples)){
dir <- Samples[i]
Total_file <- list.files(dir, recursive = T)


Noraml_file= paste(dir,Total_file[1],sep = "/")
Tumor_file <-paste(dir,Total_file[2],sep ="/")

#file1 <- "SRR7780922_Normal_Coverage_distribution.txt"
#file2 <- "SRR7780923_Tumor_Coverage_distribution.txt"

#print (Noraml_file)
#print(Tumor_file)
plotCoverage(Noraml_file, Tumor_file, 10000)

}

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






#check_bioproject[V2 >1000, ]
