info <- read.table("clipboard",sep = "\t",header = T,stringsAsFactors = F)

info[grepl("d",info$Sample.ID),]
