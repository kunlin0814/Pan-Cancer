library(tidyverse)
library(data.table)
library(gridExtra)
library(ggrepel)
library(readxl)
library(ggpubr)
library(grid)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
# gtf <- fread("G:\\MAC_Research_Data\\Pan_cancer\\Canine_Mapping_source\\Canis_familiaris.CanFam3.1.99.chr.gtf",
#              skip =5)
human_gtf <- fread("G:/MAC_Research_Data/Pan_cancer/Homo_sapiens.GRCh38.103.gtf/Homo_sapiens.GRCh38.103.gtf",
                   skip = 5)


human <- human_gtf[V3 =="gene",.(V1,V3,V4,V5,V9)]

sep_gene_stp1 <- as.data.table(str_split_fixed(human$V9," ",2))

sep_gene_stp2 <- as.data.table(str_split_fixed(sep_gene_stp1$V2,";",5))

gene_name_table <- as.data.table(str_split_fixed(sep_gene_stp2$V3,"gene_name",2))
pure_ensembl <- gsub('[\"]','',sep_gene_stp2$V1)
pure_gene <- gsub('[\ "]','',gene_name_table$V2)

human_gene_ensembl <- data.table(ensembl = pure_ensembl,
                                    gene_name = pure_gene)
fwrite(human_gene_ensembl,
       file ='G:/MAC_Research_Data/Pan_cancer/Pan_Cancer_paper/Human/OM/human_ensembl_gene.txt',
       col.names = T,row.names = F, quote = F, sep = "\t",
       eol = "\n")


a <- unique(a)
a$diff <- a$V5-a$V4
sum(a$diff)


gtf <- gtf[V1=='chr17',]

b <- gtf[grepl("NRAS",gtf$V9)]
c <- data.frame(str_split_fixed(b$V9,"gene_name",8)) 
b$V10 = b$V5-b$V4

sum(d$V10)

b$V9
a <- data.frame(str_split_fixed(gtf$V9,";",8))
b <- data.frame(str_split_fixed(a$X6,"gene_name",2)) 
b$X2
b$V9
