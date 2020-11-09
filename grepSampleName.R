library(data.table)

old <- fread("C:\\Users\\abc73_000\\Desktop\\PRJNA552034-HSA-Normal-Tumor-pair")
new <- fread("C:\\Users\\abc73_000\\Desktop\\New_PRJNA552034.txt")

normal <- new[grepl("N",new$Isolate),c("Run","Sample Name")]
tumor <- new[grepl("T",new$Isolate),c("Run","Sample Name")]

total <- new$`Sample Name`
total <- sort(total)
sample <- strsplit(total,'_')
sample_name <- NULL

for (i in 1:length(sample)){
  str <- paste(sample[[i]][1],sample[[i]][2],sep="_")
  sample_name <- c(sample_name,str)
}
sample_name <- unique(sample_name)

normal <-NULL
tumor <- NULL
sample <- NULL
for (i in total){
  if(grepl("N",i)){
    normal <- c(normal,new[`Sample Name`==i,c("Run")][,Run])
  }  
  if(grepl("T",i)){
    tumor <- c(tumor,new[`Sample Name`==i,c("Run")][,Run])
  }
  
}

combine <- cbind(sample_name,normal,tumor)

write.table(combine,file= "C:\\Users\\abc73_000\\Desktop\\new_HSA_pairs",sep="\t",
            col.names = F,
            row.names = F,
            quote = F)
    