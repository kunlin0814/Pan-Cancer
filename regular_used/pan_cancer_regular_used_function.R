library(data.table)
library(tidyverse)
library(readxl)

extract_info_by_case <- function(x){
  # use the case ID of your file to match the Case ID you want to look 
  
  index <- match(x,match_file$SampleName,nomatch = 0)
  if (index !=0){
    data <- match_file[index, .(SelfMatch,DiffFromBest)]}
  else{
    data <- data.frame(SelfMatch = NA,DiffFromBest= NA)
    
  }
  
  return(data)
}



identify_case <- function(SRR,case_total=case_total_file){
  ## this will use SRR number to match the case id
  ## need total case file to look at back
  total_row <- nrow(case_total)
  for (i in 1:total_row){
    if (SRR %in% case_total[i,]){
      case_name <- case_total[i,c("SampleName")]
      return (case_name) 
    }
    
  }
  return("NA")
  
}