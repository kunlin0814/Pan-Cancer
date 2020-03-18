library(dplyr)
library('ggplot2')
library(readxl)
library(gridExtra)



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# RNA-cat <- c('Normal', 'Tumor')
# Df  <- "G:/Pan_cancer/Pan-Cancer-Uniq_mapped_Distribution.pdf"
# for (i in RNA-cat){
# Mammary_RNA_seq <- read_excel('G:/Pan_cancer/Pan_cancer_mapping_result/CanineMC-RNAseq-Mapping.xlsx',
#                               sheet = i)
# }
# 
# ggplot(nMammary, aes(x=nMammary$Uniquely_mapped_rate)) + geom_histogram(bins = 200)+
#   xlim(min(nMammary$Uniquely_mapped_rate),max(nMammary$Uniquely_mapped_rate))+
#   xlab("Uniq_mapped_rate") +
#   ggtitle("Mammary Cancer WES")+
#   theme(axis.text=element_text(size=14),
#       axis.title=element_text(size=14,face="bold"),
#       plot.title = element_text(color = "black", size = 20, face = "bold"))
p <- list()
for(i in 1:4){
  p[[i]] <- qplot(1:10,10:1,main=i)
}
do.call(grid.arrange,p)


cat <- c('Mammary_WES_normal', 'Mammary_WES_tumor','Melanoma_WES_normal','Melanoma_WES_tumor',
         'Melanoma_WGS_normal', 'Melanoma_WGS_tumor','Osteo_WES_normal','Osteo_WES_tumor',
         'Lymphoma_WES_normal','Lymphoma_WES_tumor','Others_WES_normal','Others_WES_tumor',
         'Mammary_RNA_normal', 'Mammary_RNA_tumor')
p <- list()
#Df  <- "/Users/kun-linho/Pan-Cancer-Uniq_mapped_Distribution.pdf"
#pdf(file=Df, w=7, h=5)
for (i in cat){
  png(paste("Distribution_of_Uniquely_Mapping",i,".png",sep=""),width=6000,height=4000,res=500)
  target <- read_excel('/Volumes/Research_Data/Pan_cancer/Mamm_Melanoma_Mapping.xlsx', sheet = i)
  plot_figure <- ggplot(target, aes(x= target$Uniquely_mapped_rate)) + 
    geom_histogram(bins = 200)+
  xlim(c(min(target$Uniquely_mapped_rate)-(min(target$Uniquely_mapped_rate)*0.1), 
         max(target$Uniquely_mapped_rate)))+
  xlab("Uniq_mapped_rate")+
  ylab(" The number of the sample")+
  geom_vline(xintercept = mean(target$Uniquely_mapped_rate), linetype="dotted", 
               color = "blue", size=1.5)+
  ggtitle(i)+
  theme(axis.text=element_text(size=14),
         axis.title=element_text(size=14,face="bold"),
         plot.title = element_text(color = "black", size = 20, face = "bold"))
  p[[i]] <- plot_figure
  plot(plot_figure)
  dev.off()
 
  
}
dev.off()
#hist(x=target$Uniquely_mapped_rate, breaks = 200, 
#     xlim = c(min(target$Uniquely_mapped_rate)-0.05 ,max(target$Uniquely_mapped_rate)),
#    main = i,xlab = ("Uniq_mapped_rate"))
#Df  <- "G:/Pan_cancer/Pan-Cancer-read-pairs_Distribution.pdf"
#pdf(file=Df, w=7, h=5)
for (i in cat){
  png(paste("Distribution_of_total_pairs",i,".png",sep=""),width=6000,height=4000,res=500)
  target <- read_excel('/Volumes/Research_Data/Pan_cancer/Mamm_Melanoma_Mapping.xlsx', sheet = i)
  plot_figure <- ggplot(target, aes(x=target$Total_pairs/1000000)) + 
    geom_histogram(bins = 200)+
    #hist(x=target$Uniquely_mapped_rate, breaks = 200, 
    #     xlim = c(min(target$Uniquely_mapped_rate)-0.05 ,max(target$Uniquely_mapped_rate)),
    #    main = i,xlab = ("Uniq_mapped_rate"))
    xlim(c(min(target$Total_pairs)/1000000-min((target$Total_pairs)/1000000)*0.1,
           max(target$Total_pairs)/1000000))+
    xlab("Total Read Pairs (Million)")+
    ylab(" The number of the sample")+
    geom_vline(xintercept = mean(target$Total_pairs)/1000000, linetype="dotted", 
               color = "blue", size=1.5)+
    ggtitle(i)+
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=14,face="bold"),
          plot.title = element_text(color = "black", size = 20, face = "bold"))
  
  plot(plot_figure)
  dev.off()
  
  
}

dev.off()






