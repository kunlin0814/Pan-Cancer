library(data.table)
library(tidyverse)
library(readxl)

base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table"
seperator <- "/"
source("C:/Users/abc73/Documents/GitHub/Breed_prediction/build_sample_meta_data.R")

Path <- function(path = "clipboard") {
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

Path()


## function match_table

match_table <- function(case_id, column, table){
  info <- unlist(table[Case_ID == case_id, column, with = F])
  if (length(info)==0){
    print(case_id)
  }
  return (info)
}

## bioproject


meta <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Pan-Cancer-meta/HSA/New_PRJNA552034.txt")

path()
#### breed table 
# from breed_prediction metadata
GLM_WGS <- readClipboard()
WGS_breed <-fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Pan-Cancer-meta/GLMPRJNA579792_total.txt") 

WGS_breed_info <- WGS_breed[`Assay Type`=="WGS",.(Run,Isolate,Breed)]
GLM_WES <- readClipboard()


overlap_WGS <- intersect(GLM_WGS,GLM_WES)

WES_overlap_st <- GLM_WES %in%overlap_WGS
WGS_overlap_st <- data.frame(GLM_WGS %in%overlap_WGS)

write.table(WGS_overlap_st, file = "C:/Users/abc73/Desktop/overlap_WGS_info.txt",
            sep = "\n", col.names = T, row.names = F, quote = F)

## from excel WESQCdata
whole_table <- read.table("clipboard",sep = "\t",header = T, stringsAsFactors = F)
whole_table <- setDT(whole_table)
length(unique(whole_table$Case_ID))

pass <- unique(whole_table[The_reason_to_exclude=="Pass QC" & Case_ID!="No-Pair",])
pass <- setDT(pass)
pass_breed <- unique(pass[, .(Case_ID,Breeds)])

a <- data.frame(table(pass_breed$Breeds))
write.table(a, file = "C:/Users/abc73/Desktop/Breeds_merge_WGS_info.txt",
            sep = "\t", col.names = T, row.names = F, quote = F)


symbol_sum <- NULL
for (sample in GLM_WES){
  symbol_sum <- c(symbol_sum,unlist(whole_table[Case_ID==sample, .(Symbol)])[1])
}

write.table(symbol_sum,file = "C:/Users/abc73/Desktop/breed_WES_info.txt",
            sep = "\n",quote = F, col.names = T, row.names = F)

target <- whole_table[,.(Sample_ID,Case_ID,Tumor_Type,Breeds,Status,The_reason_to_exclude,Symbol)]
fwrite(target, file = "C:/Users/abc73/Desktop/target_WGS_info.txt",
       sep ="\t", quote = F, col.names = T, row.names = F)

### Breed prediction and assignment ##
library(readxl)
library(data.table)
# from WES_WGS_merge table
whole_wes_wgs_table <- read.table("clipboard",sep = "\t", header = T,stringsAsFactors = F)
whole_wes_wgs_table <- setDT(whole_wes_wgs_table)


breed_result <- read.table("clipboard",sep ="\t", header = T, stringsAsFactors = F)
breed_result <- setDT(breed_result)
total_samples <- unique(breed_result$SampleName)
# info <- NULL
# 
# for (i in total_samples){
#   original_label <- breed_result[SampleName==i,.(Breed)]
#   cluster <- breed_result[SampleName == i, .(BreedCluster)]
#   if (is.na(original_label)){
#     qc <- "NA"
#   }
#   else if (original_label == cluster){
#     qc <- "PassBreed"
#   }
#   else{
#     qc <- "FailBreed"
#   }
#   info <- c(info,qc)
#   
# }
# 
# breed_result$BreedQC <- info

breed_result <- assign_final_breeds(breed_result)



symbol <- NULL

for (i in total_samples){
  each <- unlist(whole_wes_wgs_table[Case_ID==i,.(Symbol)])[1]
  symbol <- c(symbol,each)
}
breed_result$Symbol <- symbol

fwrite(breed_result, file = "C:/Users/abc73/Desktop/wes_9breeds_main.txt",
       col.names = T, row.names = F, quote = F, sep = "\t",na = "NA")

assign_final_breeds <- function(meta_data){
  if ("BreedCluster" %in% colnames(meta_data)){
    total_samples <- unique(meta_data$SampleName)
    info <- NULL
    
    for (i in total_samples){
      original_label <- meta_data[SampleName==i,.(Breed)]
      cluster <- meta_data[SampleName == i, .(BreedCluster)]
      if (is.na(original_label)){
        qc <- "NA"
      }
      else if (original_label == cluster){
        qc <- "PassBreed"
      }
      else{
        qc <- "FailBreed"
      }
      info <- c(info,qc)
      
    }
    meta_data$BreedQC <- info
    final_breed <- meta_data[, "Breed"]; 
    predicted_breed_indices <- intersect(which(is.na(meta_data[, "Breed"])), which(meta_data[, "BreedCluster"] != "Unknown"));
    final_breed[predicted_breed_indices] <- meta_data[predicted_breed_indices, "BreedCluster"];
    final_breed[which(meta_data[, "BreedQC"] == "FailBreed")] <- NA;
    meta_data[, "FinalBreed"] <- final_breed;
    return (meta_data)
  }
  else {
    print("not BreedCluster assigned")
  }
}

## Breed prediction and assignment end ###

## from excel WES_BreedQCresults

total_breed_info <- read.table("clipboard",sep = "\t",header = T,stringsAsFactors = F)
total_breed_info <- setDT(total_breed_info)
WES_total_breed_info <- total_breed_info[!SampleName %in% pure_WGS, ]

write.table(WES_total_breed_info, file = "C:\\Users\\abc73\\Desktop\\wes_breed.txt",sep ="\t",
            col.names = T, row.names = F,quote = F)

# break down
target_breed <- "Unknown"
breed <- total_breed_info[Original_Breed_label ==target_breed,]
total_sample <- breed$SampleName
QC_summ <- c()
for (sample in total_sample){
  
  info <- unlist(match_table(sample,"The_reason_to_exclude",whole_qc_status))[1]
  QC_summ <- c(QC_summ, info)
}
names(QC_summ) <- NULL
breed$exclude_reason <- QC_summ
#### Original breed label
nrow(breed)
#break down
breed[,.N,keyby = .(DiseaseAcronym)]

#### add by
add_by <- total_breed_info[Original_Breed_label!=target_breed & BreedCluster==target_breed & FinalBreed==target_breed]
nrow(add_by)
add_by[,.N,keyby = .(DiseaseAcronym)]

#### remove by

remove_by <- total_breed_info[Original_Breed_label==target_breed & BreedCluster!=target_breed,]
nrow(remove_by)
# break down
remove_by[,.N,keyby = .(DiseaseAcronym)]

##### final breed label
nrow(total_breed_info[FinalBreed==target_breed,])

# break down
total_breed_info[FinalBreed==target_breed,.N, keyby=.(DiseaseAcronym)]
############## breed summary end ####

## WES
Tgen <- "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Pan-Cancer-meta/Meta/Valid_PRJNA525883_TGEN_OSA.txt"


seperator <- "/"

base_dir <- "/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table"
  #"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table"
pass_breed <- breed_info[The_reason_to_exclude =="Pass QC",]

merge_breed <- unique(pass_breed[,.(SampleName,Breed)])
merge_breed <- as.data.frame(table(merge_breed$Breed))

write.table(merge_breed, file = paste(base_dir,"merge_total_breed.txt", sep = "/"),
            sep ="\t",col.names = T,row.names = F,quote = F )

val_breed <- unique(pass_breed[Dataset=="Validation",.(SampleName,Breed)])
val_breed <- as.data.frame(table(val_breed$Breed))
sum(val_breed$Freq)
write.table(val_breed, file = paste(base_dir,"val_total_breed.txt", sep = "/"),
            sep ="\t",col.names = T,row.names = F,quote = F )
####

whole_wes_table <- fread("/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/whole_wes_table.txt")
pass_wes_table <- unique(whole_wes_table[Case_ID!="No-Pair" & The_reason_to_exclude=="Pass QC",.(Case_ID,Breed_info)])


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

fina_pass_wes_breed <- as.data.frame(table(pass_wes_table$Breed_info))
write.table(fina_pass_wes_breed, file = paste(base_dir,"final_pass_total_breed.txt", sep = "/"),
            sep ="\t",col.names = T,row.names = F,quote = F )
match_file <- fread(paste(base_dir,"tumor_normal_pair_QC_results.txt",sep = seperator))
whole_table <- fread(paste(base_dir,"whole_table.txt",sep=seperator))


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
# 
# excldue_function <- function(run){
#   if (total_file[Sample_ID == run,.(Total_pairs)<5000000]){
#     return(paste(run,"total read pairs < 5M", sep =" "))
#   }
#   else if (total_file[Sample_ID == run,.(Uniquely_coordinatly_mapped_rate)<0.6]){
#     return(paste(run,"Unique mapped_rate < 0.6"))
#   }
#   
#   else if (total_file[Sample_ID == run,.(Mean_Coverage)<30]){
#     return(paste(run,"mean coverage < 30"))
#   }
#   
#   else {
#     return("Pass QC")
#   }
# }
