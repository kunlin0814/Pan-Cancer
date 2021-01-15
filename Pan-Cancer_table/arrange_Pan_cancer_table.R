library(data.table)
library(tidyverse)
library(readxl)

excldue_function <- function(run){
  if (total_file[Sample_ID == run,.(Total_pairs)<5000000]){
    return(paste(run,"total read pairs < 5M", sep =" "))
  }
  else if (total_file[Sample_ID == run,.(Uniquely_coordinatly_mapped_rate)<0.6]){
    return(paste(run,"Unique mapped_rate < 0.6"))
  }
  
  else if (total_file[Sample_ID == run,.(Mean_Coverage)<30]){
    return(paste(run,"mean coverage < 30"))
  }
 
  else {
    return("Pass QC")
  }
}

pathPrep <- function(path = "clipboard") {
  y <- if (path == "clipboard") {
    readClipboard()
  } else {
    cat("Please enter the path:\n\n")
    readline()
  }
  x <- chartr("\\", "/", y)
  writeClipboard(x)
  return(x)
}

WP()
pathPrep()

## bioproject


meta <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Pan-Cancer-meta/HSA/New_PRJNA552034.txt")


## breed

breed <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Burair_pan_scripts/breed_prediction_test/Pan-Cancer-Breed_prediction/merge_dis_val/assignment_clusters.txt")
samples <- breed$SampleName
info <- NULL

for (i in samples){
  value <- total_file[Case_ID == i, .(Breed_info)][1]
  info <- rbind(info, value)
}

write.table(info, file = paste(base_dir,"assign_sep_breed_info.txt",sep =seperator), col.names = T, row.names = F,
            sep ="\t", quote = F)


## WGS
validation_set <- c("MT SNU", "OSA TGen", "OM Sanger", "OSA Sanger","HSA Upenn")
seperator <- "/"

base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table"

total_file <- fread(paste(base_dir,"whole_wgs_table.txt",sep = seperator))

info <- NULL 

for (i in target$V1){
  data <- total_file[Sample_ID==i, .(Sample_ID,Case_ID,Breeds,The_reason_to_exclude,Status)]
  info <- rbind(info,data)
}
write.table(info, file = paste(base_dir,"Breed_prediction_meta_WGS.txt",sep = seperator),col.names = T,
            row.names = F, quote = F, sep = "\t")

## check Breed prediction meta

meta <- read_excel("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Burair_pan_scripts/breed_prediction_test/Pan-Cancer-Breed_prediction/breed_prediction.xlsx",
                   sheet = "meta")
meta <- unique(meta[,c("SampleName","Breed")])
a <- as.data.frame(table(meta$Breed))
write.table(a, file = "C:/Users/abc73/Desktop/breed_WGS_summary.txt",col.names = T,
            row.names = F, quote = F, sep = "\t")


valid_breed_number <- as.data.frame(table(valid_breed$Breed_info))
write.table(valid_breed_number, file = paste(base_dir,"Validset_wes_breed.txt",sep = seperator),col.names = T,
            row.names = F, quote = F, sep = "\t")


total_breed <- unique(total_file[,.(Case_ID, Breed_info)])

breed <- as.data.frame(table(total_breed$Breed_info))
write.table(breed, file = paste(base_dir,"total_wes_breed.txt",sep = seperator),col.names = T,
            row.names = F, quote = F, sep = "\t")
### Breed prediction ##

library(readxl)
base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan-Cancer-Manuscript/Figure1"
breed_result <- read_excel(paste(base_dir,"WES_TableS1_1-05-21.xlsx",sep = seperator),
                           sheet = "BreedQCresults",skip = 2)


info <- NULL
  
for (i in sample_order){
  value <- total_file[Sample_ID == i, .(Symbol)]
  info <- rbind(info, value)
  }


# write.table(info, file = paste(base_dir,"QC_info.txt",sep =seperator), col.names = T, row.names = F,
#             sep ="\t", quote = F)

pass <- unique(total_file[The_reason_to_exclude=="Pass QC", .(Breed_info,Case_ID, Symbol)])

combine_breed <- as.data.frame(table(pass$Breed_info))
validation <- pass[Symbol %in% validation_set,]
validation_breed <- as.data.frame(table(validation$Breed_info))


fwrite(validation_breed, file = paste(base_dir,"validatiob_breeds.txt",sep=seperator),col.names = T,
       row.names = F, quote = F, sep = "\t")


#summary <- sapply(total_file$Sample_ID,FUN =excldue_function )
#names(summary) <- NULL

total_file$The_reason_to_exclude <- summary
output <- total_file[,.(Case_ID,Sample_ID,The_reason_to_exclude,Sample_status,Symbol)]


fwrite(output, file = paste(base_dir,"sample_status.txt",sep=seperator),col.names = T,
       row.names = F, quote = F, sep = "\t")

## WES
seperator <- "/"

base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table"

total_file <- fread(paste(base_dir,"whole_table.txt",sep = seperator))

a <- unique(total_file[Total_pairs>=5000000 & Uniquely_coordinatly_mapped_rate>=0.6 & Target_CDS_Mapping_Rates>=0.3 & Mean_Coverage<30, .(The_reason_to_exclude,Case_ID,Mean_Coverage,Symbol)])
a$Case_ID
                
#&Uniquely_coordinatly_mapped_rate>=0.6 & Target_CDS_Mapping_Rates >=0.3 & Mean_Coverage>=30 ] 
a <- a[!Callable_bases %in% c("Normal_sample","No-Pair")]
a$Callable_bases <- as.numeric(a$Callable_bases)
b <- a[a$Callable_bases<10000000]

a$Callable_bases <- as.numeric(a$Callable_bases)

validation <- unique(total_file[,.(Case_ID,Breed_info)])
breed <- as.data.frame(table(validation$Breed_info)) 
fwrite(breed, file= paste(base_dir,"all_breed.txt",sep=seperator),col.names = T, row.names = F,
       quote = F, sep = "\t")  

unique(total_file$Symbol)

tgen_sample <- fread(paste(base_dir,"tgne_WGS.txt",sep = seperator),header = F)
tgen_meta <- fread(paste(base_dir,"TGENbreedsOSA.txt",sep=seperator))

breed <- tgen_meta[`Sample Name` %in% tgen_sample$V1,.(Breed)]
fwrite(breed, file= paste(base_dir,"Tgne_WGS_breed.txt",sep = seperator),col.names = T,row.names = F,
       quote = F, sep="\t")


match_file <- fread(paste(base_dir,"tumor_normal_pair_QC_results.txt",sep = seperator))
whole_table <- fread(paste(base_dir,"whole_GLM_tgne_cases.txt",sep=seperator),header = F)
breed_info <- fread(paste(base_dir,"pan-cancer_breeds.txt",sep=seperator))
breed_info[,Original_label_breeds]
sample_status <- fread(paste(base_dir,"Pan-cancer_sample.txt",sep=seperator),header = F)

extract_info_by_case <- function(x){
  index <- match(x,match_file$SampleName,nomatch = 0)
  if (index !=0){
  data <- match_file[index, .(SelfMatch,DiffFromBest)]}
  else{
    data <- data.frame(SelfMatch = NA,DiffFromBest= NA)
    
  }
  
  return(data)
}


extract_breed_info_by_case <- function(x){
  index <- match(x,breed_info[,Case_ID],nomatch = 0)
  if (index !=0){
    data <- breed_info[index, .(Original_label_breeds)]}
  else{
    data <-NA
    
  }
  
  return(data)
}

extract_sample_status_info_by_case <- function(x){
  index <- match(x,sample_status[,V1],nomatch = 0)
  if (index !=0){
    data <- sample_status[index, .(V2)]}
  else{
    data <-"NA"
    
  }
  
  return(data)
}



#brd_info <- unlist(sapply(whole_table$Case_ID,FUN = extract_breed_info_by_case))

sample_status_info <- unlist(sapply(whole_table$V1,FUN = extract_sample_status_info_by_case))
whole_table$sample_status <- sample_status_info

  fwrite(whole_table, file = paste(base_dir,"whole_sample_status.txt",sep=seperator),quote = F,
         col.names = T, row.names = F, sep = "\t", na = "NA")



whole_table$Breed_info <- brd_info
whole_table$Sample_status <- sample_status_info
fwrite(whole_table, file = paste(base_dir,"whole_new_table.txt",sep=seperator),quote = F,
       col.names = T, row.names = F, sep = "\t", na = "NA")

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

