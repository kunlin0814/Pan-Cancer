library(data.table)


source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
out_put_base <- "C:/Users/abc73/Desktop"
meta <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Pan-Cancer-meta/Melanoma_PRJEB12081_Meta.txt")

final_out <- meta[,.(Run,BioSample,title)]

normal_sample <- final_out[grepl('b$',final_out$title),]
tumor_sample <- final_out[grepl('a$',final_out$title),]
meta_sample <- final_out[grepl('[c$]',final_out$title),]
meta2_sample <- final_out[grepl('[d$]',final_out$title),]
normal_sample <- normal_sample[order(title,Run)]
tumor_sample <- tumor_sample[order(title,Run)]
meta_sample <- meta_sample[order(title,Run)]

final_normal <- NULL
final_tumor <- NULL

for (index in 1:length(normal_sample$title)){
  each_normal <- normal_sample$title[index]
  each_tumor <- tumor_sample$title[index]
  normal_summ <- normal_sample[title ==each_normal,.SD,keyby = .(title)]
  each_normal_table <- dcast(normal_summ, title~Run,value.var = 'Run')
  colnames(each_normal_table) <- c("sample_names","run1","run2")
  
  tumor_summ <- tumor_sample[title ==each_tumor,.SD,keyby = .(title)]
  each_tumor_table <- dcast(tumor_summ, title~Run,value.var = 'Run')
  colnames(each_tumor_table) <- c("sample_names","run1","run2")
  final_normal <- rbindlist(list(final_normal,each_normal_table))
  final_tumor <- rbindlist(list(final_tumor,each_tumor_table))
}
final_normal <- unique(final_normal)
final_tumor <- unique(final_tumor)


final_meta <- NULL

for (index in 1:length(meta_sample$title)){
  each_meta <- meta_sample$title[index]
  meta_summ <- meta_sample[title ==each_meta,.SD,keyby = .(title)]
  each_meta_table <- dcast(meta_summ, title~Run,value.var = 'Run')
  colnames(each_meta_table) <- c("sample_names","run1","run2")
  final_meta <- rbindlist(list(final_meta,each_meta_table))
}

final_meta <- unique(final_meta)

fwrite(final_normal,paste(out_put_base,"OM_normal_sample.txt",sep = "/"),
       col.names = T, row.names = F, quote = F, sep = "\t",
       eol = "\n")

fwrite(final_tumor,paste(out_put_base,"OM_tumor_sample.txt",sep = "/"),
       col.names = T, row.names = F, quote = F, sep = "\t",
       eol = "\n")
fwrite(final_meta,paste(out_put_base,"OM_meta_sample.txt",sep = "/"),
       col.names = T, row.names = F, quote = F, sep = "\t",
       eol = "\n")


fwrite(meta2_sample,paste(out_put_base,"OM_meta2_sample.txt",sep = "/"),
       col.names = T, row.names = F, quote = F, sep = "\t",
       eol = "\n")
