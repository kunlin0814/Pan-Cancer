### This script analyzed the mutation among two subtypes of MC and use wilcox test to identify 
### signifncat mutation betweens two subtypes

library(dplyr)
library(readxl)
library(Biobase)
library(FSA)
library(data.table)
library(ggplot2)

SimpleAdenoma <- read_excel("C:\\Users\\abc73_000\\Desktop\\MC_VAF\\Total_MC_VAF.xlsx",
                            sheet = "Simple_Ade")

SimpleCarcionma <-read_excel("C:\\Users\\abc73_000\\Desktop\\MC_VAF\\Total_MC_VAF.xlsx",
                             sheet = "Simple_Car")


# Mutation <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Total_MC_VAF.xlsx",
#                        sheet = "Total_Mut")


SimpleAdenoma <- as.matrix(SimpleAdenoma)
SimpleCarcionma <- as.matrix(SimpleCarcionma)
Overlap_Genes  <- colnames(SimpleAdenoma)
colNumber <- ncol(SimpleAdenoma)
both_empty_column <- c()

for (i in 2:colNumber) {
  SimpleAdenoma_column_value <- as.numeric(SimpleAdenoma[, i])
  SimpleCarcionma_column_value <- as.numeric(SimpleCarcionma[, i])
  if ( all(as.numeric(SimpleAdenoma[, i])==0) & all(as.numeric(SimpleCarcionma[, i]) ==0)) {
    both_empty_column <- c(both_empty_column, i)
  }
}

if (!is.null(both_empty_column)){
  SimpleAdenoma <- SimpleAdenoma[,-both_empty_column]
  SimpleCarcionma <- SimpleCarcionma[,-both_empty_column]
  } else {
  SimpleAdenoma <- as.matrix(SimpleAdenoma)
  SimpleCarcionma <- as.matrix(SimpleCarcionma)
  }

Overlap_Genes  <- colnames(SimpleAdenoma)
#Overlap_Genes <- Overlap_Genes[-1]
SimpleAdenoma_list <- list()
SimpleCarcionma_list <- list()
#Overlap_Genes1 <- colnames(SimpleCarcionma)
for (i in 1:length(Overlap_Genes)){
  SimpleAdenoma_list[[Overlap_Genes[i]]] = SimpleAdenoma[,i]
  SimpleCarcionma_list[[Overlap_Genes[i]]] = SimpleCarcionma[,i]
}

## to get the plot for the species distribution and for each species distribution
Df  <- paste("distribution", "_","Gene_location","_","MC_Mutation.pdf", sep="", collapse="");
pdf(file=Df, w=7, h=5)
par( mar=c(2.1,4.1,2.1,1.1) )
layout(m=matrix(1:2, 2, 1))

SimpleCarcionma_significant_total <- list()
SimpleAdenoma_significant_total <- list()
significant_Gene_name <- c()
non_sig_Gene_name <- c()
SimpleAdenoma_sign_value <- c()
SimpleCarcionma_sign_value <- c()
p_value_Gene <- c()
fisher_pvalue_Gene <- c()
for (i in 2:length(Overlap_Genes)){
  SimpleAdenoma_Gene_value <- as.numeric(SimpleAdenoma[, i])
  SimpleCarcionma_Gene_value <- as.numeric(SimpleCarcionma[, i])
  Gt5SampleSimpleAdenoma <- length(which(SimpleAdenoma_Gene_value>0))
  Gt5SampleSimpleCarcionma <- length(which(SimpleCarcionma_Gene_value>0))
  wilcox <- wilcox.test(SimpleAdenoma_Gene_value,SimpleCarcionma_Gene_value, alt="two.sided",paired = F, correct=T)
  Count0SimpleAdenoma <- length(SimpleAdenoma_Gene_value)-Gt5SampleSimpleAdenoma
  Count0SimpleCarcionma <- length(SimpleCarcionma_Gene_value)-Gt5SampleSimpleCarcionma
  fisher_matrix <- matrix(c(Gt5SampleSimpleAdenoma,Gt5SampleSimpleCarcionma,Count0SimpleAdenoma,Count0SimpleCarcionma),nrow=2, ncol=2)
  fisher_test <- fisher.test(fisher_matrix)
  fisher_p_value <- fisher_test$p.value
  p_value=wilcox$p.value
  if (Gt5SampleSimpleAdenoma >=5 |Gt5SampleSimpleCarcionma >=5){ # here the p_value is the raw pvalue of wilcox text
    significant_Gene <- Overlap_Genes[i]
    p_value_Gene <- c(p_value_Gene, p_value)
    #print(fisher_p_value)
    #print(sprintf("Gt5SampleSimpleAdenoma is %s%s",  Gt5SampleSimpleAdenoma,significant_Gene))
    #print(sprintf("Gt5SampleSimpleCarcionma is %s", Gt5SampleSimpleCarcionma ))
    SimpleAdenoma_sign_value <- as.numeric(SimpleAdenoma[, i])
    SimpleCarcionma_sign_value <- as.numeric(SimpleCarcionma[, i])
    SimpleAdenomaframe <- data.frame(value = SimpleAdenoma_sign_value, Status = c("SimpleAdenoma") )
    SimpleCarcionmaframe <- data.frame(value = SimpleCarcionma_sign_value, Status = c("SimpleCarcinoma"))
    MCdata <- rbind(SimpleAdenomaframe,SimpleCarcionmaframe)
    #,SimpleCarcinoma= SimpleCarcionma_sign_value )
    SimpleAdenoma_significant_total[[significant_Gene]] <- SimpleAdenoma_sign_value
    p <- ggplot(data = MCdata, aes(x=factor(Status,levels = c("SimpleAdenoma", "SimpleCarcinoma")),
               y=as.numeric(value),fill=Status,color=Status)) +
      
      geom_point(size=1.6,shape=20,position = position_jitterdodge(jitter.width = 0.28)) +
      ylab("VAF")+
      ggtitle(paste("MC" ,paste(significant_Gene,sep="_", collapse="")))
    print(p)
    #hist_value <- SimpleAdenoma_significant_total[[significant_Gene]]
    #hist(as.numeric(hist_value),breaks = 100, xlab = 'Value',main = paste("SimpleAdenoma of" ,paste(significant_Gene,sep="_", collapse="")))

    SimpleCarcionma_significant_total[[significant_Gene]] <- SimpleCarcionma_sign_value
    #hist_value1 <-SimpleCarcionma_significant_total[[significant_Gene]]
    #hist(as.numeric(hist_value1), breaks = 100,xlab = 'Enrichemnt',main = paste("SimpleCarcionma of" ,paste(significant_Gene,sep="_", collapse="")))
    fisher_pvalue_Gene <- c(fisher_pvalue_Gene,fisher_p_value )
    significant_Gene_name <- c(significant_Gene_name, Overlap_Genes[i] )
  }else{
    non_sig_Gene_name <- c(non_sig_Gene_name,Overlap_Genes[i])
    }
}

dev.off()




Gene_pvalue <- data.frame(wilcoxPvalue=p_value_Gene, fisherPvalue=fisher_pvalue_Gene, Gene= significant_Gene_name )



#write.table(Gene_pvalue, file ="/Users/kun-linho/Desktop/Gene_pvalue.txt", row.names=F, sep ="\t",quote = F  )

### Once we have the statistically significant Gene, we need to adjust the pvalue ####
Data = Gene_pvalue[order(Gene_pvalue$wilcoxPvalue),]
headtail(Data)
adjust_pvalue <- p.adjust(Data$wilcoxPvalue, method ="hochberg", n = length(Data$Gene))
Data$wilcoxBH = p.adjust(Data$wilcoxPvalue, method = "BH",n = length(Data$Gene))
Data$fisherBH = p.adjust(Data$fisherPvalue, method = "BH",n = length(Data$Gene))
#Data$Bonferroni = p.adjust(Data$wilcoxPvalue, method = "bonferroni")
# Data$Holm = p.adjust(Data$ wilcoxPvalue, method = "holm")
# Data$Hochberg = p.adjust(Data$ wilcoxPvalue, method = "hochberg")
# Data$Hommel = p.adjust(Data$ wilcoxPvalue, method = "hommel")
# Data$BY = p.adjust(Data$ wilcoxPvalue, method = "BY")


#write.table(Data, file ="/Users/kun-linho/Desktop/Gene_adjust_different_methods.txt", row.names=F, sep ="\t",quote = F  )
data_with_adjust_pvalue <- data.frame(wilcoxBH_adjust_p= Data$wilcoxBH, 
                                      fisherBH_adjust_p= Data$fisherBH,
                                      Gene=Data$Gene, 
                                      wilcoxRow_P= Data$wilcoxPvalue,
                                      fisherRow_p= Data$fisherPvalue)
### data_with_adjust_pvalue contains row p value and adjust p value ####
data_with_adjust_pvalue <- data_with_adjust_pvalue[order(data_with_adjust_pvalue$Gene), ]
write.table(data_with_adjust_pvalue,
            file ="C:\\Users\\abc73_000\\Desktop\\Pure_SignificantMutation_MC.txt", row.names=F, sep ="\t",quote = F)



wilcox.test(as.numeric(SimpleAdenoma[,"PIK3CA"]),
                       as.numeric(SimpleCarcionma[,"PIK3CA"]),
                                  alt="two.side",paired = F)


median(as.numeric(SimpleAdenoma[,"PIK3CA"]))
median(as.numeric(SimpleCarcionma[,"PIK3CA"]))


SimpleAdenoma_significant_Gene_list <- list()
SimpleCarcionma_significant_Gene_list <- list()

for (i in Gene_pvalue$Gene){
  SimpleAdenoma_significant_Gene_list[[i]]= SimpleAdenoma_list[[i]]
  SimpleCarcionma_significant_Gene_list[[i]]= SimpleCarcionma_list[[i]]
}

SimpleAdenoma_list_dataframe <-as.data.frame(SimpleAdenoma_significant_Gene_list)
SimpleCarcionma_list_dataframe <- as.data.frame(SimpleCarcionma_significant_Gene_list)
#data.frame(matrix(unlist(SimpleAdenoma_significant_Gene_list), nrow = max(lengths(SimpleAdenoma_significant_Gene_list)),byrow=T),stringsAsFactors=F)
#SimpleCarcionma_list_dataframe <- data.frame(matrix(unlist(SimpleCarcionma_significant_Gene_list), nrow = max(lengths(SimpleCarcionma_significant_Gene_list)), byrow=T),stringsAsFactors=F)
#colnames(SimpleAdenoma_list_dataframe) <- Gene_pvalue$Gene
#colnames(SimpleCarcionma_list_dataframe) <- Gene_pvalue$Gene

## here we need significant Gene read as input for each sample ##
#SimpleCarcionma_sig_data <- read_excel("/Users/kun-linho/Desktop/Gene_SimpleAdenoma_SimpleCarcionma_distribution.xlsx", sheet = 'SimpleCarcionma_Gene_significant')
#SimpleAdenoma_sig_data <- read_excel("/Users/kun-linho/Desktop/Gene_SimpleAdenoma_SimpleCarcionma_distribution.xlsx", sheet = 'SimpleAdenoma_Gene_significant')
col_names <-colnames(SimpleAdenoma_list_dataframe)
SimpleCarcionma_sig_data <- as.matrix(SimpleCarcionma_list_dataframe)
SimpleAdenoma_sig_data <- as.matrix(SimpleAdenoma_list_dataframe)
#colnames(SimpleCarcionma_sig_data) <- NULL
#colnames(SimpleAdenoma_sig_data) <- NULL
SimpleCarcionma_sample_amount <- length(SimpleCarcionma_sig_data[,1])
SimpleAdenoma_sample_amount <- length(SimpleAdenoma_sig_data[,1])

