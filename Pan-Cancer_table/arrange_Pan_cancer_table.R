library(data.table)
library(tidyverse)
library(readxl)

## WES
seperator <- "/"

base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table"

match_file <- fread(paste(base_dir,"tumor_normal_pair_QC_results.txt",sep = seperator))
whole_table <- fread(paste(base_dir,"whole_table.txt",sep=seperator))


extract_info_by_case <- function(x){
  index <- match(x,match_file$SampleName,nomatch = 0)
  if (index !=0){
  data <- match_file[index, .(SelfMatch,DiffFromBest)]}
  else{
    data <- data.frame(SelfMatch = NA,DiffFromBest= NA)
    
  }
  
  return(data)
}


pair_info <- as.data.frame(sapply(whole_table$Case_ID,FUN = extract_info_by_case))
pair_info <- as.data.frame(t(pair_info))
fwrite(pair_info,file = paste(base_dir,"total_pair_info.txt",sep =seperator), col.names = T,
       row.names = T,
       sep ="\t",quote = F,na = "NA")


## WGS ## file up case names
seperator <- "/"
base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/WGS_analysis"
case_total_file <- read_excel(paste(base_dir,"WGS_QC.xlsx",sep=seperator),sheet = "Total_cases")
table_file <- read_excel(paste(base_dir,"WGS_QC.xlsx",sep=seperator),sheet = "Sheet1")

identify_case <- function(SRR,case_total=case_total_file){
  total_row <- nrow(case_total)
  for (i in 1:total_row){
  if (SRR %in% case_total[i,]){
    case_name <- case_total[i,c("SampleName")]
    return (case_name) 
  }
  
  }
  return("NA")

}

case_name <- unlist(sapply(table_file$file_name,identify_case))
table_file$Case_ID <- case_name

write.table(table_file,file = paste(base_dir,"WGS_final_table.txt",sep=seperator),
            sep = "\t",quote = F,col.names = T,row.names = F)

