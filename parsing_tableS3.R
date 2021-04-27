library(readxl)
library(data.table)
library(janitor)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
setwd('G:/MAC_Research_Data/need_to_revised/FigureS3')

tableS3 <- read_excel("TableS3_4-22-21.xlsx",
                      skip = 1, sheet = "Mutation+CNA")

tableS3 <- data.frame(tableS3, stringsAsFactors = F)
first_part <- tableS3[1:3,]
first_part <- as.data.frame(t(first_part),stringsAsFactors = F)
setDT(first_part, keep.rownames = TRUE,check.names=T)
names(first_part) <- unlist(as.character(first_part[1,]))
first_part <- first_part[-1,]
second_part <- tableS3[18:nrow(tableS3),]
second_part <- as.data.frame(t(second_part),stringsAsFactors = F)
setDT(second_part, keep.rownames = TRUE)
names(second_part) <- unlist(as.character(second_part[1,]))
second_part <- second_part[-1,]
total_table <- left_join(first_part, second_part,by= "Sample")
sum(!is.na(total_table[,.(PIK3CA)]))
#total_table[is.na(total_table)] <- 0


which(tableS3$Sample == "PIK3CA")


a = tableS3[c(1,9360),which(!is.na(tableS3[9360,]))]
data <- as.data.frame(t(a))
data <- setDT(data)

colnames(data) <- c("Disease","PIK3CA")

b= data[grepl("H1047R",data$PIK3CA) & Disease=="MT",]

data <- fread("/Volumes/Research/MAC_Research_Data/need_to_revised/Manuscript/TableS3_4-20-21_MT_H1047R.csv")
a = data[Disease=="MT",.(PIK3CA)]
b = a[grepl("H1047R",a$PIK3CA),]
