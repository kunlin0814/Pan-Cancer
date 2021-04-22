library(readxl)

tableS3 <- read_excel("/Volumes/Research/MAC_Research_Data/need_to_revised/FigureS3/TableS3_4-20-21.xlsx",
                      skip = 1, sheet = "Mutation+CNA")
tableS3 <- setDF(tableS3)

tableS3[,1]

which(tableS3$Sample == "PIK3CA")
a = tableS3[c(1,9360),which(!is.na(tableS3[9360,]))]
data <- as.data.frame(t(a))
data <- setDT(data)

colnames(data) <- c("Disease","PIK3CA")

b= data[grepl("H1047R",data$PIK3CA) & Disease=="MT",]

data <- fread("/Volumes/Research/MAC_Research_Data/need_to_revised/Manuscript/TableS3_4-20-21_MT_H1047R.csv")
a = data[Disease=="MT",.(PIK3CA)]
b = a[grepl("H1047R",a$PIK3CA),]
