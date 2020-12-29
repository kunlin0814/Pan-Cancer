library(tidyverse)
library(readxl)
#library(wesanderson)
library(RColorBrewer)
library(data.table)

## the script contains function that use cases name to identify the subtype of lymphoma dataset in pan-cancer study

total_file <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\Var_Lofreq_Somatic_seq\\Somatic_seq.xlsx",
                         sheet ='PON_DBSNP_prioritize')
LYMcellSubtype <- fread("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\source\\LymSubtype.txt")
LYMcellSubtype$SubType

exclude <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Figure1\\Original_Data_summary.xlsx",
                      sheet ="Total_excluded")

check_case <- function(SRR){
  total_row <- nrow(val_sample)
  for (i in 1:total_row){
    if (SRR %in% val_sample[i, ]){
      
      return (val_sample[i,c("Cases") ]) }
    
  }
  return ("NaN")
}


check_Lymsubtype <- function(SampleName){
  value <- match(SampleName, LYMcellSubtype$SampleName,nomatch=0)
  if (value!=0){
    return (LYMcellSubtype[value,c("SubType") ])
  }
  else{
    return ("NaN")
  }
}
subtype <- unlist(sapply(total_file$file_name,check_Lymsubtype))
total_file <- cbind(total_file, subtype) %>% 
  write.table("C:\\Users\\abc73_000\\Desktop\\WithBcellsubtype.txt",
              sep ='\t',row.names = F,quote = F)




origin_breeds <- sapply(total_file$Sample_ID,check_CaseID)