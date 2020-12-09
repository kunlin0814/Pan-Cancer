library(tidyverse)
library(data.table)
library(gridExtra)
library(ggrepel)
library(readxl)
library(ggpubr)
library(grid)

excldue <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Original_Data_summary.xlsx",
                      sheet = "Before_Matching_excluded")

regular.text <- element_text(colour="black",size=20)

my_jitter <- function(dataframe, dot_size=1, abs_text_size=16, xangle=30) {
  regular.text <- element_text(colour="black",size=abs_text_size);
  xaxis_just <- ifelse(xangle > 0, 1, 0.5);
  
  p <- ggplot(dataframe, aes(x=x, y=y)) + geom_jitter(color="black", size=dot_size, shape=20, position=position_jitter(width=0.2));
  p <- p + stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", size=0.2, width=0.8, color = "red", show.legend=FALSE);
  p <- p + theme(axis.text=regular.text, axis.title=regular.text, axis.text.x = element_text(angle=xangle, hjust=xaxis_just), legend.position="none",
                 panel.background=element_blank(), axis.line=element_line(color="black"),
                 plot.title=element_text(face="plain", size=abs_text_size, hjust=0.5));
  return(p);
}

clean_table <- function(table){
  ref_snp <- which(nchar(table$ref)==1)
  alt_snp <- which(nchar(table$alt)==1)
  keep_index <- intersect(ref_snp,alt_snp)
  final_table <- table[keep_index,]
  final_table$mut_type <- paste(final_table$ref,final_table$alt,sep =">")
  
  return(final_table)
}

convert_mutation_type <- function(x) {
  if(x == "A>C") {
    return("T>G");
  } else if(x == "A>G") {
    return("T>C");
  } else if(x == "A>T") {
    return("T>A");
  } else if(x == "G>A") {
    return("C>T");
  } else if(x == "G>C") {
    return("C>G");
  } else if(x == "G>T") {
    return("C>A");
  } else {
    return(x);
  }
}

count_mutation_types <- function(clean_data, signature_levels) {
  result <- rep(0, length(signature_levels));
  names(result) <- signature_levels;
  total_mut <- sum(clean_data$number)
  df <- setNames(data.frame(matrix(ncol = 6, nrow = 1)), signature_levels)
  #df <- setDT(df)
  if (total_mut ==0){
    
    final_result <- result
  }
  else{
    for ( i in signature_levels){
      if (!i %in% clean_data$conver_mut_type){
        df[, i] <-  0  
      }
      else{
        df[, i] <- sum(clean_data[conver_mut_type==i,.(number)])
      }
    }
    for ( i in signature_levels){
      result[i] <- df[,i]
    }
    final_result <- result/total_mut
  }
  return(final_result);
}


my_barplot <- function(dataframe, fill_colors, abs_text_size=16, xangle=30) {
  regular.text <- element_text(colour="black",size=abs_text_size);
  xaxis_just <- ifelse(xangle > 0, 1, 0.5);
  
  p <- ggplot(dataframe, aes(x=x, y=y, fill=fill)) + geom_bar(stat="identity", position="stack") + 
    scale_fill_manual(values=fill_colors);
  p <- p + theme(axis.text=regular.text, axis.title=regular.text, axis.text.x = element_text(angle=xangle, hjust=xaxis_just),
                 panel.background=element_blank(), axis.line=element_line(color="black"),
                 legend.title=regular.text, legend.position="top", legend.text=regular.text,
                 plot.title=element_text(face="plain", size=abs_text_size, hjust=0.5));
  return(p);
}

check_empty_fill <- function(data_table){
  col_name <- c("number","ref","alt")
  if (nrow(data_table)==0){
    data_table <- setNames(data.frame(matrix(ncol = 3, nrow = 1)), col_name)  
    data_table$number <- 0
    data_table$ref <- "C"
    data_table$alt <- "A"
  }
  else{
    data_table <- data_table
    colnames(data_table) <- col_name
  }
  return(data_table)
}


prepare_bar_plot_from_table <- function(y,sample_name){
  signatures_fill <- factor(names(y), signature_levels);
  x <- rep(sample_name, length(signature_levels));
  data <- data.frame(x=x,y=y,fill=signatures_fill)
  
  return(data)
}


signature_colors <- c("cyan","black","red","gray","green","pink");
signature_levels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");


cancer_type <- c("MC", "OM")
attribute <- c("before", "after")

base <- "C:\\Users\\abc73_000\\Desktop\\bases_sub"

Cancer_type <- "OM"
#C:\\Users\\abc73_000\\Desktop\\bases_sub\\MC\\Mutect1\\CMT-2\\before\\CMT-2_Mutect1_before5steps.txt
#C:\\Users\\abc73_000\\Desktop\\bases_sub\\MC\\Mutect1\\CMT-2\\after\\CMT-2_Mutect1_after5steps.txt
OM_sample <-  list.files (path = "C:\\Users\\abc73_000\\Desktop\\bases_sub\\OM\\Mutect1")
MC_sample <- list.files (path = "C:\\Users\\abc73_000\\Desktop\\bases_sub\\MC\\Mutect1")
col_name <- c("number","ref","alt")
xangle <- 45


pdf("C:\\Users\\abc73_000\\Desktop\\bases_sub\\OM_data_6bases.pdf"
    , height=8.98, width=8.84);

for (sample in OM_sample){
  M1before_file_name <- paste(sample,"Mutect1_before5steps.txt",sep = '_')
  M1after_file_name <- paste(sample,"Mutect1_after5steps.txt",sep = '_')
  M2before_file_name <- paste(sample,"Mutect2_before5steps.txt",sep = '_')
  M2after_file_name <- paste(sample,"Mutect2_after5steps.txt",sep = '_')
  
  Mutect1_before <-fread(paste(base,Cancer_type,"Mutect1",sample,"before",M1before_file_name ,sep="\\"))
  Mutect1_after <- fread(paste(base,Cancer_type,"Mutect1",sample, "after",M1after_file_name ,sep="\\"))
  Mutect2_before <-fread(paste(base,Cancer_type,"Mutect2",sample,"before",M2before_file_name ,sep="\\"))
  Mutect2_after <- fread(paste(base,Cancer_type,"Mutect2",sample, "after",M2after_file_name ,sep="\\"))
  
  
  Mutect1_before <- check_empty_fill(Mutect1_before)
  Mutect1_after <- check_empty_fill(Mutect1_after)
  Mutect2_before <- check_empty_fill(Mutect2_before)
  Mutect2_after <- check_empty_fill(Mutect2_after)
  
  
  clean_Mutect1_before <- clean_table(Mutect1_before)
  clean_Mutect1_after <- clean_table(Mutect1_after)
  clean_Mutect2_before <- clean_table(Mutect2_before)
  clean_Mutect2_after <- clean_table(Mutect2_after)
  
  clean_Mutect1_before$conver_mut_type <- sapply(clean_Mutect1_before$mut_type, convert_mutation_type)
  clean_Mutect1_after$conver_mut_type <- sapply(clean_Mutect1_after$mut_type, convert_mutation_type)
  clean_Mutect2_before$conver_mut_type <- sapply(clean_Mutect2_before$mut_type, convert_mutation_type)
  clean_Mutect2_after$conver_mut_type <- sapply(clean_Mutect2_after$mut_type, convert_mutation_type)
  
  clean_Mutect1_before <- count_mutation_types(clean_Mutect1_before, signature_levels)
  clean_Mutect1_after <- count_mutation_types(clean_Mutect1_after, signature_levels)
  clean_Mutect2_before <- count_mutation_types(clean_Mutect2_before, signature_levels)
  clean_Mutect2_after <- count_mutation_types(clean_Mutect2_after, signature_levels)
  
  
  
  data1 <- prepare_bar_plot_from_table(clean_Mutect1_before,sample)
  data2 <- prepare_bar_plot_from_table(clean_Mutect1_after,sample)
  data3 <- prepare_bar_plot_from_table(clean_Mutect2_before,sample)
  data4 <- prepare_bar_plot_from_table(clean_Mutect2_after,sample)
  #  
  
  p1 <- my_barplot(data1, fill_colors=signature_colors, xangle=xangle);
  p1 <- p1 + labs(title="Mutect1_Before", y="", x="", fill="Mut.");
  
  p2 <- my_barplot(data2, fill_colors=signature_colors, xangle=xangle);
  p2 <- p2 + labs(title="Mutect1_after", y="", x="", fill="Mut.");
  
  p3 <- my_barplot(data3, fill_colors=signature_colors, xangle=xangle);
  p3 <- p3 + labs(title="Mutect2_Before", y="", x="", fill="Mut.");
  
  p4 <- my_barplot(data4, fill_colors=signature_colors, xangle=xangle);
  p4 <- p4 + labs(title="Mutect2_after", y="", x="", fill="Mut.");
  #
  # 
  # ggarrange(p1, p2, p3 + p4,
  #            labels = c("A", "B", "C","D"),
  #            ncol = 2, nrow = 2)
  grid.arrange(p1, p2,p3, p4,nrow = 2, top = textGrob(sample,gp=gpar(fontsize=24,font=3)))
  
}
dev.off()


my_barplot(data1,fill_colors=signature_colors, abs_text_size=16, xangle=30)

plot(ggarrange(plotlist=append(dotplots, barplots), nrow=2, ncol=4))





x <- rep("test", length(signature_levels));
data <- data.frame(x=x,y=y,fill=signatures_fill)

xangle <- 45
p <- my_barplot(data, fill_colors=signature_colors, xangle=xangle);
p <- p + labs(title="", y="", x="", fill="Mut.");
barplots[["apple"]] <- p;
print(p)

