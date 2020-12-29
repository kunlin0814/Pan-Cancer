library(ggforce)
library(tidyr) 
library(dplyr)
library(ggplot2)
library(ggforce)
fill_colors = signature_colors
x <- mut_before_count_sum$x
y <- mut_before_count_sum$y;
mutation_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");
fill <- mut_before_count_sum$fill
fill <- factor(fill, levels=mutation_types);
samples <- unique(x);
sample_order <- order(sapply(samples, function(s) {sum(y[which(x == s)])}), decreasing=T);
x <- factor(x, levels=samples[sample_order]);
data <- data.frame(x=as.numeric(x), y=y, fill=fill);

pdf(paste(OM_base,"test_same_order_OM_samples_data_6bases_count.pdf",sep="\\")
    , height=12.94, width=12.94);

ggplot(data, aes(x = as.numeric(x), y = y, fill = fill))+
  geom_bar(stat = "identity", position = "dodge", width = 0.7)+
  scale_fill_manual(values=fill_colors) + 
  scale_x_continuous(
    breaks = 1:length(levels(data$x)),
    label = levels(data$x)
  )+
  facet_zoom(x= x %in% c("DD0001","DD0002"), y<300,horizontal=F)+
  labs(x = "", y = " ")+
  theme(legend.background = element_rect(fill="transparent"),
        legend.title = element_blank(), axis.text.x=element_text(angle=55, vjust=1,  hjust=1,size = 10))

dev.off()




