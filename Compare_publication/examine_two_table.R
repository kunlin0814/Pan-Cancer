library(data.table)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
source("C:/Users/abc73/Documents/GitHub/VAF/six_base_function_util.R")
seperator <- "/"
base_dir <- "C:/Users/abc73/Desktop/re-run_OM"
new <- fread(paste(base_dir,"Sanger_OM_pair_spread.txt",sep =seperator ))
old <- fread(paste(base_dir,"combined_melanoma.txt",sep =seperator ),header = F)
old$sample_names <- sapply(old$V1,convert_sample)


new$combine_normal <- paste(new$normal1,new$normal2,sep = "-")
new$combine_tumor <- paste(new$tumor1,new$tumor2,sep = "-")

new_target <- new[,.(`#sample_names`,combine_normal,combine_tumor)]
old_target <- old[,.(V1,V2,V3)]
colnames(new_target)[1] <- "sample_names"
colnames(old_target) <- c("sample_names","combine_normal","combine_tumor")
setkey(old_target,sample_names)
old_target[new_target,nomatch=0]

a = setdiff(old_target,new_target)

fwrite(old_target, file = paste(base_dir,"old_OM_complete_pair.txt",sep = "/"),
       col.names = T,row.names = F,quote = F, sep = "\t")


new_target[old_target==new_target]
