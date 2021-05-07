##### Pan cancer Supp Fig. 6b
library(ggplot2)
library(MASS)
library(scales)
library(data.table)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
setwd("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate")
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
whole_wes_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt") 
exclude <- unique(unlist(whole_wes_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

mutect2 <- fread("Before5steps_total_TMB_Mutect2.txt")
mutect2 <- mutect2[!file_name %in% exclude, ]
mutect2$subtype <- match_vector_table(mutect2$file_name, "DiseaseAcronym2", whole_wes_table)
mutect2 <-mutect2[subtype!="UCL"]

atLeastTwo <- fread("at_least_two_no_UCL_excludeQC_TMB_04_15.txt")
mutect2$callable <- match_info_table(mutect2$file_name,"sample_names", "callable",atLeastTwo, string_value = F)

mutect2$TMb <- (mutect2$NumberofMut*1000000)/mutect2$callable

fwrite(mutect2, file= "mutect2_tmb.txt",
       col.names = T, row.names = F, sep = "\t", eol = "\n",quote = F)
#mutect2 = read.delim('Mutect2_Burair_Bias3_no_UCL_excludeQC_TMB_03_21.txt',header = T)
#atLeastTwo = read.delim('at_least_two_no_UCL_excludeQC_TMB_03_21.txt',header = T)


mutect2$tumor_type = factor(mutect2$subtype,levels = c('MT','GLM','BCL','TCL','OM','OSA','HSA'))
mutect2$TMb = as.numeric(as.vector(mutect2$TMb))
abs_text_size <- 12
xangle <- 45
regular.text <- element_text(size=abs_text_size);
png(file = "FigS6B_mutect2.png", width = 1200, height = 1000, units = "px", res = 500)

p <- ggplot(mutect2, aes(x=tumor_type, y=log10(TMb+0.01))) + 
  geom_jitter(aes(colour = tumor_type), shape=16, position=position_jitterdodge(),size = 0.5) + 
  ylab("TMB")+
  ggtitle("MuTect2")+
  theme_classic()+ 
  scale_color_manual(values=c("darkblue","darkblue","darkblue","firebrick3","firebrick3","firebrick3","firebrick3")) +
  #scale_y_continuous(trans='log10',breaks = c(0,0.1,1,10))+
  coord_cartesian(ylim = c(-2,1))+
  theme(axis.text=regular.text, 
        axis.title=regular.text, 
        axis.text.x = element_text(angle=xangle, hjust=xaxis_just,colour = c("darkblue","darkblue","darkblue","firebrick3","firebrick3","firebrick3","firebrick3")),
        axis.text.y = element_text(colour = 'black'),
        axis.title.x =element_blank(),
        axis.title.y= element_text(colour = 'black'),
        legend.position="none")+
  stat_summary(geom = "crossbar", width=0.7, fatten=2, color="black", fun.data = function(x){c(y=mean(x), ymin=mean(x), ymax=mean(x))},position = position_dodge(0.78))
  
print(p)



dev.off()



atLeastTwo$tumor_type = factor(atLeastTwo$tumor_type,levels = c('MT','GLM','BCL','TCL','OM','OSA','HSA'))
atLeastTwo$TMb = as.numeric(as.vector(atLeastTwo$TMb))


png(file = "FigS6B_atLeastTwo.png", width = 1200, height = 1000, units = "px", res = 500)

ggplot(atLeastTwo, aes(x=tumor_type, y=log10(TMb+0.01))) + 
  geom_jitter(aes(colour = tumor_type), shape=16, position=position_jitterdodge(),size = 0.5) + 
  ylab("TMB")+
  ggtitle("At least two software")+
  theme_classic()+ 
  coord_cartesian(ylim = c(-2,1))+
  #scale_y_continuous(trans='log10',breaks = c(0,0.1,1,10))+
  scale_color_manual(values=c("darkblue","darkblue","darkblue","firebrick3","firebrick3","firebrick3","firebrick3")) +
  stat_summary(geom = "crossbar", width=0.7, fatten=2, color="black", fun.data = function(x){c(y=mean(x), ymin=mean(x), ymax=mean(x))},position = position_dodge(0.78))+
  theme(axis.text=regular.text, 
        axis.title=regular.text, 
        axis.text.x = element_text(angle=xangle, hjust=xaxis_just,colour = c("darkblue","darkblue","darkblue","firebrick3","firebrick3","firebrick3","firebrick3")),
        axis.text.y = element_text(colour = 'black')##### Pan cancer Supp Fig. 6b
        
library(ggplot2)
library(MASS)
library(scales)
library(data.table)
source("C:/Users/abc73/Documents/GitHub/R_util/my_util.R")
setwd("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Mutation_rate_VAF/Mut_rate")

whole_wes_table <- fread("G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/all_pan_cancer_wes_metatable_04_09.txt") 
exclude <- unique(unlist(whole_wes_table[The_reason_to_exclude!="Pass QC",.(Case_ID)]))

mutect2 <- fread("Before5steps_total_TMB_Mutect2.txt")
mutect2 <- mutect2[!file_name %in% exclude, ]
mutect2$subtype <- match_vector_table(mutect2$file_name, "DiseaseAcronym2", whole_wes_table)
mutect2 <-mutect2[subtype!="UCL"]

atLeastTwo <- fread("at_least_two_no_UCL_excludeQC_TMB_04_15.txt")
mutect2$callable <- match_info_table(mutect2$file_name,"sample_names", "callable",atLeastTwo, string_value = F)

mutect2$TMb <- (mutect2$NumberofMut*1000000)/mutect2$callable
#mutect2 = read.delim('Mutect2_Burair_Bias3_no_UCL_excludeQC_TMB_03_21.txt',header = T)
#atLeastTwo = read.delim('at_least_two_no_UCL_excludeQC_TMB_03_21.txt',header = T)


mutect2$tumor_type = factor(mutect2$subtype,levels = c('MT','GLM','BCL','TCL','OM','OSA','HSA'))
mutect2$TMb = as.numeric(as.vector(mutect2$TMb))
abs_text_size <- 12
xangle <- 45
regular.text <- element_text(size=abs_text_size);
png(file = "FigS6B_mutect2.png", width = 1200, height = 1000, units = "px", res = 500)

p <- ggplot(mutect2, aes(x=tumor_type, y=log10(TMb+0.01))) + 
  geom_jitter(aes(colour = tumor_type), shape=16, position=position_jitterdodge(),size = 0.5) + 
  ylab("TMB")+
  ggtitle("Mutect2")+
  theme_classic()+ 
  scale_color_manual(values=c("darkblue","darkblue","darkblue","firebrick3","firebrick3","firebrick3","firebrick3")) +
  #scale_y_continuous(trans='log10',breaks = c(0,0.1,1,10))+
  coord_cartesian(ylim = c(-2,1))+
  theme(axis.text=regular.text, 
        axis.title=regular.text, 
        axis.text.x = element_text(angle=xangle, hjust=xaxis_just,colour = c("darkblue","darkblue","darkblue","firebrick3","firebrick3","firebrick3","firebrick3")),
        axis.text.y = element_text(colour = 'black'),
        axis.title.x =element_blank(),
        axis.title.y= element_text(colour = 'black'),
        legend.position="none")+
  stat_summary(geom = "crossbar", width=0.7, fatten=2, color="black", fun.data = function(x){c(y=mean(x), ymin=mean(x), ymax=mean(x))},position = position_dodge(0.78))
  
print(p)



dev.off()



atLeastTwo$tumor_type = factor(atLeastTwo$tumor_type,levels = c('MT','GLM','BCL','TCL','OM','OSA','HSA'))
atLeastTwo$TMb = as.numeric(as.vector(atLeastTwo$TMb))


png(file = "FigS6B_atLeastTwo.png", width = 1200, height = 1000, units = "px", res = 500)

ggplot(atLeastTwo, aes(x=tumor_type, y=log10(TMb+0.01))) + 
  geom_jitter(aes(colour = tumor_type), shape=16, position=position_jitterdodge(),size = 0.5) + 
  ylab("TMB")+
  ggtitle("At least two software")+
  theme_classic()+ 
  coord_cartesian(ylim = c(-2,1))+
  #scale_y_continuous(trans='log10',breaks = c(0,0.1,1,10))+
  scale_color_manual(values=c("darkblue","darkblue","darkblue","firebrick3","firebrick3","firebrick3","firebrick3")) +
  stat_summary(geom = "crossbar", width=0.7, fatten=2, color="black", fun.data = function(x){c(y=mean(x), ymin=mean(x), ymax=mean(x))},position = position_dodge(0.78))+
  theme(axis.text=regular.text, 
        axis.title=regular.text, 
        axis.text.x = element_text(angle=xangle, hjust=xaxis_just,colour = c("darkblue","darkblue","darkblue","firebrick3","firebrick3","firebrick3","firebrick3")),
        axis.text.y = element_text(colour = 'black'),
        axis.title.x =element_blank(),
        axis.title.y= element_text(colour = 'black'),
        legend.position="none")


dev.off()




