library(tidyverse)
library(data.table)
library(gridExtra)
library(ggrepel)
library(readxl)
library(ggpubr)
library(grid)


regular.text <- element_text(colour="black",size=20)

convert_MT_sample <- function(MT_sample){
  each <- str_split(MT_sample,"-")
  while (nchar(each[[1]][2])<3){
    each[[1]][2] <- paste("0",each[[1]][2], sep="")
  }
  return (paste(each[[1]][1],each[[1]][2],sep = "-"))
}



convert_sample <- function(sample){
  if(grepl("-1",sample)){
    split_words <- str_split(sample,'-')[[1]][1]
    sanger_sample <- paste(split_words,'c',sep = "")
  }else if(grepl("-2",sample)){
    split_words <- str_split(sample,'-')[[1]][1]
    sanger_sample <- paste(split_words,'d',sep = "")
  }
  else{
    sanger_sample <- paste(sample,'a',sep = "")
  }
  return (sanger_sample)
  
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

conver_to_six_bases_basedon_mutation <- function(data_frame, sample,cancer_type,rates = T){
  fontsize <- 20;
  signature_colors <- c("cyan","black","red","gray","green","pink");
  signature_levels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
  
  data_frame$mut_type <- paste(data_frame$ref,data_frame$alt,sep = ">")
  data_frame$conver_mut_type <- sapply(data_frame$mut_type, convert_mutation_type)
  data_frame <- setDT(data_frame)
  data_info <- as.data.frame(table(data_frame[sample_name == sample,.(conver_mut_type)]))
  colnames(data_info) <- c("conver_mut_type","number")
  data_info <- setDT(data_info)
  if (rates){
    data_info <- count_mutation_rates(data_info,signature_levels)
    data_info <- prepare_bar_plot_from_table(data_info,sample)
    data_info$tumor_type <- "OM"
    rownames(data_info) <- NULL
  }else{
    data_info <- count_mutation_number(data_info,signature_levels)
    data_info <- prepare_bar_plot_from_table(data_info,sample)
    data_info$tumor_type <- "OM"
    rownames(data_info) <- NULL  
  }  
  
  return (data_info)
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
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size=fontsize,face="plain",colour="black"),
      legend.title= element_blank(), legend.text = element_text(size=fontsize,face="plain",colour="black"));
  return(p)
}
