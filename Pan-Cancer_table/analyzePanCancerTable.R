library(tidyverse)
library(readxl)
library(wesanderson)
library(RColorBrewer)
library(data.table)

check_coverage <- function(x){
  if (x < 30){
    return ("<30")
  } else if (x >30 & x<50){
    return ("30<x<50")}
  else if(x >50 & x <100){
     return ("50<x<100")
   }else if(x>=100){
     return("x>=100")
   }
}

total_data <- read_excel("C:\\Users\\abc73_000\\Desktop\\New_WES_QC_dataset.xlsx")
total_data <- setDT(total_data)
stat <- sapply(total_data$mean,check_coverage) 
total_data$coverage_stat <- stat
sample <- total_data$Case_ID

tumor_stat <- NULL
for (samp in sample){
  each <- total_data[Case_ID == samp, ]
  each_stat <- each[Status=="Tumor",.(coverage_stat)]
  tumor_stat <- c(tumor_stat,each_stat)
}
tumor_stat <- unlist(tumor_stat)

total_data$tumor_coverage <- tumor_stat

unique(total_data$coverage_stat)


both_stat <- NULL


# 
# each <- total_data[Case_ID == "004",]
# each_coverage <- unlist(each[,.(coverage_stat)])

for (samp in sample){
  each <- total_data[Case_ID == samp,]
  each_coverage <- each$coverage_stat

if (each_coverage[1] == each_coverage[2]){
  status = each_coverage[1]
}else{
  if("<30" %in% each_coverage){
    status <- "<30"
  }
  else if("30<x<50" %in% each_coverage){
    status <- "30<x<50"
  }
  else if ("50<x<100" %in% each_coverage){
    status <- "50<x<100"
  }
  else{
    status <- "x>=100"
  }
  
}
  both_stat <- c(both_stat,status)
}
total_data$both_status <- both_stat

total_tumor <- total_data[Status =="Tumor",]


sample_dist <- total_data[,.(number = length(Sample_ID)),keyby = .(Tumor_Type,Status)]
lt30 <- total_data[mean <30,.(number = length(Sample_ID)),keyby = .(Tumor_Type,Status)]  
sum(lt30$number)


lt50 <- total_data[mean <50,.(number = length(Sample_ID)),keyby = .(Tumor_Type,Status)]
sum(lt50$number)

gt50 <- total_data[mean >=50,.(number = length(Sample_ID)),keyby = .(Tumor_Type,Status)]  


gt50_100 <- total_data[mean >=50 & mean<100,.(number = length(Sample_ID)),keyby = .(Tumor_Type,Status)] 
sum(gt50_100$number)
gt100 <- total_data[mean >=100,.(number = length(Sample_ID)),keyby = .(Tumor_Type,Status)] 
sum(gt100$number)
