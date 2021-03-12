library(data.table)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")


mt <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Pan-Cancer-meta/Meta/Valid_PRJNA552905_MC.txt")
mt <- mt[,.(Run,Age)]


osa <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Pan-Cancer-meta/Meta/Valid_PRJNA525883_TGEN_OSA.txt")
osa <- osa[,.(Run,Age)]

hsa <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Pan-Cancer-meta/Meta/Valid_PRJNA417727_HSA.txt")
hsa <- hsa[,.(Run,Age)]


total <- rbindlist(list(mt, osa, hsa))

total_breed_info <- read.table("clipboard",sep = "\t",header = T,stringsAsFactors = F)

#target_id, table_column_to_use, output_column ,table, string_value = T

age <- match_info_table(total_breed_info$Sample_ID,"Run", "Age",total)

total_breed_info$Age <- age

fwrite(total_breed_info, file = "C:\\Users\\abc73\\Desktop\\validation_data_set_age.txt",
       col.names = T, row.names = F, quote = F, sep = "\t")
