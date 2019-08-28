# Import the packages
library(Biobase)
library(NMF)

# Read in -  NMF work
file=read.table(file='/scratch/kh31516/Mammary/Mammary_cancer_hallmark.txt',header=T,sep="\t",stringsAsFactors = FALSE)
# Setting the first column as row names
row.names(file)=as.character(file[,1])
file=file[,-1]
#class(file)
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
#write.table(coln,"coln.txt",quote=F,row.names=F,col.names=F)
# Re-read in the annotation file after been annotated
#name=read.table(file='/Users/kun-linho/Desktop/coln_treat_filtered copy.txt',header=T)
#summary(name)
#all(rownames(name)==colnames(file_nonneg)) #check samples name
#metadata=data.frame(labelDescription=c('Status','Location'),row.names=c("Status","Location"))
phenoData=new("AnnotatedDataFrame", data=coln)
#              data=coln, varMetadata=metadata)
#phenoData


# experiment description
experimentData=new("MIAME",name="Kun-Lin",title="Hallmark NMF Tumor vs Normal") 


# assembling an expression-set
exampleSet=ExpressionSet(assayData=exprs,
                         experimentData=experimentData,
                         annotation='hgu95av2')
                         #phenoData=phenoData,)



# Get the best rank based on cophenetic parameters - by running

estim.r=nmf(exampleSet,2:6,nrun=200,seed=123456,.options='P10')
tiff("NMF-ranks-hallmark_genes-Mammary_cancer.tiff", width = 4000, height = 2000, res=500)
plot(estim.r)
dev.off()

tiff("NMF-ConsensusMaps-hallmark_genes-Mammary_cancer.tiff", width = 6500, height = 4000, res=300)
consensusmap(estim.r,annCol=exampleSet,labCol=NA,labRow=NA)
dev.off()

# Run NMF on the best rank # need to re-run based on the classification
est=nmf(exampleSet,2,nrun=200,seed=123456,.options='P10')


tiff("NMF-ConsensusMaps-rank2-hallmark_genes-Mammary_cancer.tiff", width = 2000, height = 2000, res=300)
# Consensus map (rank = 3)
consensusmap(est,annCol=exampleSet,labCol=colnames(file),labRow=NA)
dev.off()


tiff("NMF-Basis-rank2-hallmark_genes-Mammary_cancer.tiff", width = 4000, height = 2000, res=300)
# basismap & coefmap
#layout(cbind(1,2))
basismap(est,subsetRow=F)
dev.off()
#a=phenoData$Status
#b=phenoData$Location
#c=data.frame(a,b)
tiff("NMF-CoefficientMaps-rank2-hallmark_genes-Mammary_cancer.tiff", width = 4000, height = 2000, res=300)
coefmap(est) # annCol can ignore
dev.off()