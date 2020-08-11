# Import the packages
#source('http://bioconductor.org/biocLite.R')
#biocLite('NMF', siteRepos=c('http://web.cbio.uct.ac.za/~renaud/CRAN'), type='both')

#library(devtools)
#install_github('renozao/pkgmaker', ref = 'develop')
#install_github('renozao/NMF',ref= 'devel')
#library(devtools)
library(Biobase)
library(NMF)
#file=read.delim('/scratch/kh31516/Mammary/source/Mammary_CRC_matrix.txt',header=T, sep='\t')
# Read in -  NMF work
file=read.table('/scratch/kh31516/Mammary/source/Mammary_Melanoma_hallmark.txt',header=T, sep='\t')
row.names(file)=as.character(t(file[,1]))
file=file[,-1]
coli0 <- which(colSums(file) == 0)
coli_na <- which(colSums(is.na(file)) > 0)
rowi0 <- which(rowSums(file)==0)
row_na <- which(rowSums(is.na(file))>0)
#file=file[c(-rowi0), -coli0]
file=as.data.frame(t(file))

#file <- as.matrix(file)
class(file)
#which(rowSums(file[,1:171])==0)
#which(rowSums(is.na(file[,1:171]))>0)

#i0 <- which(colSums(file) == 0)
#i_na <- which(colSums(is.na(file)) > 0)
#file=file[,-1]
#discard "cob.N14.77.N" "X869T1A" "X869T1B"
#file=file[,c(-15,-16,-26)]

# the values for file[10608:10618,ncol(file)] is -inf, which is under the cell line X257191
# Regard them as missing values and use interval DEA approach to fill in all these blanks
#for (i in 10608:10618){
#  file[i,ncol(file)]<-NA
#  minimun=min(file[i,],na.rm=T)
#  maximun=max(file[i,],na.rm=T)
#  if (is.na(file[i,ncol(file)])){
#    file[i,ncol(file)]=runif(1,minimun,maximun)
#  }
#} 

# Normalized the dataset as we transform the values in the input matrix to positive values by adding 
# the absolute value of the lowest value to all cells
#which(file==min(file))
file_nonneg=file+abs(min(file))
exprs=as.matrix(file_nonneg)


# prepare the phenotype dataset
coln=colnames(file_nonneg)
coln=data.frame(coln)
#write.table(coln,"/scratch/kh31516/Mammary/source/Mamm_melan_hallamrkannot.txt",quote=F,row.names=F,col.names=F)
# Re-read in the annotation file after been annotated
name=read.table(file='/scratch/kh31516/Mammary/source/Mamm_melan_hallamrkannot.txt',header=T)
summary(name)
all(rownames(name)==colnames(file_nonneg)) #check samples name
metadata=data.frame(labelDescription='Cancer_type',row.names="Cancer_type")
phenoData=new("AnnotatedDataFrame", data=name,varMetadata=metadata )
#              data=coln, varMetadata=metadata)
#phenoData


# experiment description
experimentData=new("MIAME",name="Kun-Lin",title="Mammary_Melanoma NMF Tumor vs Normal") 


# assembling an expression-set
exampleSet=ExpressionSet(assayData=exprs,
                         experimentData=experimentData,
                         annotation='hgu95av2',
                         phenoData=phenoData)
#exampleSet


# file=read.delim('/scratch/kh31516/Mammary/source/Mammary_CRC_matrix.txt',header=T, sep='\t')
# Get the best rank based on cophenetic parameters - by running
estim.r=nmf(exampleSet,2:6,nrun=500,seed=123456, .options='P10')
tiff("500 run NMF-ranks-Hallmark-Mammary_Melanoma.tiff", width = 4000, height = 2000, res=500)
plot(estim.r)
dev.off()
tiff("500 run NMF-ranks-Hallmark-ConsensusMaps-Mammary_Melanoma.tiff", width = 6500, height = 4000, res=300)
consensusmap(estim.r,annCol=exampleSet,labCol=NA,labRow=NA)
dev.off()
