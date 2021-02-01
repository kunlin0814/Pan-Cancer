library(data.table)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
gtf <- fread("C:/Users/abc73/Desktop/TestTMB/Canis_familiaris.CanFam3.1.99.chr.gtf",
             skip = 5,stringsAsFactors = F)
candidate_list <- fread("C:/Users/abc73/Desktop/TestTMB/total_target_ensmbl_id.txt", header = F)
gtf_gene <- gtf[V3=="gene",]

sep_gene_stp1 <- as.data.table(str_split_fixed(gtf_gene$V9," ",2))

sep_gene_stp2 <- as.data.table(str_split_fixed(sep_gene_stp1$V2,";",3))

pure_ensembl <- gsub('[\"]','',sep_gene_stp2$V1)

gtf_gene$ensembl <- pure_ensembl


output <- gtf_gene[,.(ensembl,V1,V4,V5)]
col_names <- c("Ensemble_id","chromosome_location",	"gene_start","gene_end")

colnames(output) <- col_names

final_out <- output[Ensemble_id %in% candidate_list$V1,][order(Ensemble_id)]

fwrite(final_out, file = "C:/Users/abc73/Desktop/TestTMB/Pan_cancer_sign_target_ensembl_location.txt",
       quote = F, col.names = T, row.names = F, sep = "\t",eol = "\n")


