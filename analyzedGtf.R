library(tidyverse)
library(data.table)
library(gridExtra)
library(ggrepel)
library(readxl)
library(ggpubr)
library(grid)

gtf <- fread("G:\\MAC_Research_Data\\Pan_cancer\\Canine_Mapping_source\\Canis_familiaris.CanFam3.1.99.chr.gtf",
             skip =5)


a <- gtf[V3 =="exon",.(V1,V3,V4,V5)]
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
