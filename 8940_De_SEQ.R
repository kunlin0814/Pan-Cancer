# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("DESeq2")

Sys.setenv(LANG = "en")
library("DESeq2")

# Set the working directory (where you downloaded your count files)
directory <- "G:/MAC_Research_Data/UGA/2Spring/8940/HW_8940_RNA_seq"
# Set the prefix for each output file name
outputPrefix <- "binf8940_DESeq2_parsed"

sampleFiles<-
  c("UHR_Rep1_gene.tsv","UHR_Rep2_gene.tsv","UHR_Rep3_gene.tsv","HBR_Rep1_gene.tsv","HBR_Rep2_gene.tsv","HBR_Rep3_gene.tsv")

sampleNames <- c("UHR 1","UHR 2","UHR 3","HBR 1","HBR 2","HBR 3")
sampleCondition <- c("UHR","UHR","UHR","HBR","HBR","HBR")
sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)
treatments = c("UHR","HBR")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design = ~ condition)

colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

# DESeq2 run
dds <- DESeq(ddsHTSeq)
res <- results(dds)
# Order results by padj value (most significant to least)
res= subset(res, padj<0.05)
res <- res[order(res$padj),]
# Save data results and normalized reads to csv
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
names(resdata)[1] <- 'gene'
write.csv(resdata, file = paste0(outputPrefix, "-results-with-normalized_parsed.csv"))
# Produce DataFrame of results of statistical tests
mcols(res, use.names = T)
# Replace outlier value with estimated value as predicted by distrubution using "trimmed mean" approach.
# Recommended if you have several replicates per treatment. DESeq2 will automatically do this 
# if you have 7 or more replicates
ddsClean <- replaceOutliersWithTrimmedMean(dds)
ddsClean <- DESeq(ddsClean)
tab <- table(initial = results(dds)$padj < 0.05,
             cleaned = results(ddsClean)$padj < 0.05)
addmargins(tab)
resClean <- results(ddsClean)
resClean = subset(res, padj<0.05)
resClean <- resClean[order(resClean$padj),]
#PLOTING (also in Rstudio)
# MA plot of RNAseq data for entire dataset with genes with padj < 0.05 colored Red
plotMA(ddsClean, ylim=c(-8,8),main = "RNAseq experiment")

# Transform raw counts into normalized values (DESeq2 has two options: 1) rlog transformed and 2) variance
# stabilization variance stabilization is very good for heatmaps!)

rld <- rlogTransformation(ddsClean, blind=T)
vsd <- varianceStabilizingTransformation(ddsClean, blind=T)
# heatmap of data
library("RColorBrewer")
library("gplots")
#library(heatmap2)
library(heatmap3)
# Top 50 expressed genes with heatmap3
plot_result <- png("C:\\Users\\abc73_000\\Desktop\\heatmap_gene_expression.png",
                   width=2800,height=1800,res=350)  
par(mar=c(1,1,1,1))
select <- order(rowMeans(counts(ddsClean,normalized=T)),decreasing=T)[1:50]

my_palette <- colorRampPalette(c("blue",'black','yellow'))(n=72)

heatmap.2(assay(vsd)[select,], col=my_palette,
         scale="row", key=T, keysize=1, symkey=T,
         density.info="none", trace="none",dendrogram = "both",
         cexCol=0.6, labRow=F,
         main="Expressed Genes Heatmap")

dev.off()


