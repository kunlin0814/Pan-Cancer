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

count_mutation_rates <- function(clean_data, signature_levels) {
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


count_mutation_number <- function(clean_data, signature_levels) {
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
    final_result <- result
  }
  return(final_result);
}



my_barplot <- function(dataframe, fill_colors, abs_text_size=12, xangle=30) {
  regular.text <- element_text(colour="black",size=abs_text_size);
  xaxis_just <- ifelse(xangle > 0, 1, 0.5);
  
  p <- ggplot(dataframe, aes(x=x, y=y, fill=fill)) + geom_bar(stat="identity", position="stack") + 
    scale_fill_manual(values=fill_colors);
  p <- p + theme(axis.text=regular.text, axis.title=regular.text, axis.text.x = element_text(angle=xangle, hjust=xaxis_just),
                 panel.background=element_blank(), axis.line=element_line(color="black"),
                 legend.title=regular.text, 
                 legend.position="top", 
                 legend.text=regular.text,
                 plot.title=element_text(face="plain", size=abs_text_size, hjust=0.5));
  return(p);
}

check_empty_fill <- function(data_table){
  
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
  
  #p <- my_barplot(data, fill_colors=signature_colors, xangle=xangle);
  return(data)
}

find_six_base <- function(data,rate = T,sample_name){
  data <- check_empty_fill(data)
  data <- clean_table(data)
  data$conver_mut_type <- sapply(data$mut_type, convert_mutation_type)
  
  if(rate==T){
    data <- count_mutation_rates(data,signature_levels)
  }
  else{
    data <- count_mutation_number(data,signature_levels)
    
  }
  data <- prepare_bar_plot_from_table(data,sample_name)
  return(data)
}

sample_signature <- function(data,fill_colors,title){
  x <- data$x
  y <- data$y;
  mutation_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
  fill <- data$fill
  fill <- factor(fill, levels=mutation_types);
  samples <- unique(x);
  sample_order <- order(sapply(samples, function(s) {sum(y[which(x == s)])}), decreasing=T);
  x <- factor(x, levels=samples[sample_order]);
  data <- data.frame(x=x, y=y, fill=fill);
  
  p <- ggplot(data, aes(x=x, y=y, fill=fill)) + geom_bar(stat="identity",position='stack', width=0.6)+
    ggtitle(title)+ 
    scale_fill_manual(values=fill_colors)+
    theme(
      plot.title = element_text(size = 20, face = "bold"),
      axis.title.x = element_blank(),
      #element_text(face="plain",colour="black",size=fontsize),
      axis.title.y = element_blank(),
      #element_text(face="plain",colour="black",size=fontsize),
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=fontsize,face="plain",colour="black"),
      legend.title= element_blank(), legend.text = element_text(size=fontsize,face="plain",colour="black"));
  return(p)
}


fontsize <- 20;
signature_colors <- c("cyan","black","red","gray","green","pink");
signature_levels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
col_name <- c("number","ref","alt")

base <- "C:\\Users\\abc73_000\\Desktop\\bases_sub"

Cancer_type <- "OM"
OM_sample <-  sort(list.files (path = "C:\\Users\\abc73_000\\Desktop\\bases_sub\\OM\\Mutect1"))
MC_sample <- sort(list.files (path = "C:\\Users\\abc73_000\\Desktop\\bases_sub\\MC\\Mutect1"))
MC_sample <- MC_sample[-match('CMT-33', MC_sample)]

col_name <- c("number","ref","alt")
xangle <- 45

# pdf("C:\\Users\\abc73_000\\Desktop\\bases_sub\\try_combined_MC_data_6bases_count.pdf"
#     , height=8.94, width=12.84);

  if (Cancer_type =="OM"){
  total_sample <-OM_sample
  }else{
    total_sample <- MC_sample
  }
mut_before_count_sum <- NULL
mut_after_count_sum <- NULL
mut2_before_count_sum <- NULL
mut2_after_count_sum <- NULL

for (sample in total_sample){
  M1before_file_name <- paste(sample,"Mutect1_before5steps.txt",sep = '_')
  M1after_file_name <- paste(sample,"Mutect1_after5steps.txt",sep = '_')
  M2before_file_name <- paste(sample,"Mutect2_before5steps.txt",sep = '_')
  M2after_file_name <- paste(sample,"Mutect2_after5steps.txt",sep = '_')
  
  Mutect1_before <-fread(paste(base,Cancer_type,"Mutect1",sample,"before",M1before_file_name ,sep="\\"))
  Mutect1_after <- fread(paste(base,Cancer_type,"Mutect1",sample, "after",M1after_file_name ,sep="\\"))
  Mutect2_before <-fread(paste(base,Cancer_type,"Mutect2",sample,"before",M2before_file_name ,sep="\\"))
  Mutect2_after <- fread(paste(base,Cancer_type,"Mutect2",sample, "after",M2after_file_name ,sep="\\"))
  
  

  # mut_before_rate <- find_six_base(Mutect1_before,sample,rate =T);
  # mut_after_rate <- find_six_base(Mutect1_after,sample,rate =T);
  # mut2_before_rate <- find_six_base(Mutect2_before,sample,rate =T);
  # mut2_after_rate <- find_six_base(Mutect2_after,sample,rate =T);

  
  mut_before_count <- find_six_base(Mutect1_before,rate =F,sample_name=sample);
  mut_after_count <- find_six_base(Mutect1_after,rate =F,sample_name=sample);
  mut2_before_count <- find_six_base(Mutect2_before,rate =F,sample_name=sample);
  mut2_after_count <- find_six_base(Mutect2_after,rate =F,sample_name=sample);
  
  mut_before_count_sum <- rbind(mut_before_count_sum,mut_before_count)
  mut_after_count_sum <- rbind(mut_after_count_sum,mut_after_count)
  mut2_before_count_sum <- rbind(mut2_before_count_sum,mut2_before_count)
  mut2_after_count_sum <- rbind(mut2_after_count_sum,mut2_after_count)
  
  
  ## Plot the signature
#   p1 <- my_barplot(mut_before_rate, fill_colors, abs_text_size=12, xangle=30)
#   p1 <- p1 + labs(title="mut_before_Ratio", y="", x="", fill="Mut.");
#   p2 <- my_barplot(mut_after_rate,fill_colors, abs_text_size=12, xangle=30)
#   p2 <- p2 + labs(title="mut_after_Ratio", y="", x="", fill="Mut.");
#   p3 <- my_barplot(mut2_before_rate,fill_colors, abs_text_size=12, xangle=30)
#   p3 <- p3 + labs(title="mut2_before_Ratio", y="", x="", fill="Mut.");
#   p4 <- my_barplot(mut2_after_rate,fill_colors, abs_text_size=12, xangle=30)
#   p4 <- p4 + labs(title="mut2_after_Ratio", y="", x="", fill="Mut.");
#   
#   
#   p1_1 <- my_barplot(mut_before_count,fill_colors, abs_text_size=12, xangle=30)
#   p1_1 <- p1_1 + labs(title="mut_before_Counts", y="", x="", fill="Mut.");
#   p2_1 <- my_barplot(mut_after_count,fill_colors, abs_text_size=12, xangle=30)
#   p2_1 <- p2_1 + labs(title="mut_after_Counts", y="", x="", fill="Mut.");
#   p3_1 <- my_barplot(mut2_before_count,fill_colors, abs_text_size=12, xangle=30)
#   p3_1 <- p3_1 + labs(title="mut2_before_Counts", y="", x="", fill="Mut.");
#   p4_1 <- my_barplot(mut2_after_count,fill_colors, abs_text_size=12, xangle=30)
#   p4_1 <- p4_1 + labs(title="mut2_after_Counts", y="", x="", fill="Mut.");
#   
#   lay <- rbind(c(1,2,3,4),
#                c(5,6,7,8))
#   #c(6,7,8,9,9))
#   grid.arrange(p1,p1_1,p2,p2_1,
#                p3,p3_1,p4,p4_1,
#                nrow = 2, top = textGrob(sample,gp=gpar(fontsize=24,font=3)),layout_matrix = lay)
# }
}

## check what samples are not normal
ex <- setDT(mut_before_count_sum)
answer <- ex[, .(total=sum(y)),keyby = .(x)][order(-total)]


# dev.off()

## All the samples together
pdf(paste("C:\\Users\\abc73_000\\Desktop\\bases_sub\\",Cancer_type,"_samples_data_6bases_count.pdf",sep="")
    , height=12.94, width=12.94);

p1 <- sample_signature(mut_before_count_sum,fill_colors = signature_colors, title = "Mutect Before")
p2 <- sample_signature(mut_after_count_sum,fill_colors = signature_colors, title = "Mutect After")
p3 <- sample_signature(mut2_before_count_sum,fill_colors = signature_colors, title = "Mutect2 Before")
p4 <- sample_signature(mut2_after_count_sum,fill_colors = signature_colors, title = "Mutec2 After")


 grid.arrange(p1,p2,
             nrow = 2, top = textGrob(Cancer_type,gp=gpar(fontsize=24,font=3)))

 grid.arrange(p3,p4,
             nrow = 2, top = textGrob(Cancer_type,gp=gpar(fontsize=24,font=3)))

dev.off()

# Sanger Data comparison
sanger_data <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\OM_mutation_compare_with_Sanger\\Sanger_mutation.xlsx",
                          sheet ="Canine", skip = 29)

sanger_data <- setDT(sanger_data)

pure_sig_data <- sanger_data[,.(Ref,Alt,Sample)]
colnames(pure_sig_data) <- c("ref","alt","sample")
clean_sanger <- clean_table(pure_sig_data)
clean_sanger$conver_mut_type <- sapply(clean_sanger$mut_type, convert_mutation_type)

a <- clean_sanger[sample=="DD0001a",]
b <- data.frame(table(a$conver_mut_type))
b <- setDT(b)
colnames(b) <- c("conver_mut_type","number")
c <- count_mutation_number(b,signature_levels )
prepare_bar_plot_from_table(c,)
find_six_base(b,rate =F,sample_name=sample);

result <- rep(0, length(signature_levels));
names(result) <- signature_levels;
df <- setNames(data.frame(matrix(ncol = 6, nrow = 1)), signature_levels)

summary <- NULL
for (samp in total_sample){
  
  a <- clean_sanger[sample==samp,]
  b <- data.frame(table(a$conver_mut_type))
  #b$sample_name <- samp
  b <- setDT(b)
  colnames(b) <- c("conver_mut_type","number")
  c <- count_mutation_number(b,signature_levels )
  d <- prepare_bar_plot_from_table(c,samp)
  summary <- rbind(summary, d)
}


pdf(paste("C:\\Users\\abc73_000\\Desktop\\bases_sub\\","Sanger","_samples_data_6bases_count.pdf",sep="")
    , height=12.94, width=12.94);




sanger <- sample_signature(summary,fill_colors = signature_colors, title = "Sanger_Data")
grid.arrange(p2,sanger,
             nrow = 2, top = textGrob(Cancer_type,gp=gpar(fontsize=24,font=3)))

grid.arrange(p4,sanger,
             nrow = 2, top = textGrob(Cancer_type,gp=gpar(fontsize=24,font=3)))

dev.off()
# sanger compaiosn end 
