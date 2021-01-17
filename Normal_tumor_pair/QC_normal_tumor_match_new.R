library("ggplot2");
library("ggthemes");
library("RColorBrewer");

cancer_types <- c("MT", "GLM", "LYM", "OM", "OSA", "HSA", "UCL");
#dataset_order <- c("MT CUK", "MT SNU", "MT UGA", "GLM JL", "LYM BI", "OM SC1", "OM SC2", "OSA BI", "OSA TGen", "OSA SC", "HSA BI", "HSA UPenn", "UCL BI");
# OM Sanger and OSA Sanger all failed sequencing QC, MT UGA was removed
dataset_order <- c("MT CUK", "MT SNU", "GLM JL", "LYM BI", "OM SC1", "OSA BI", "OSA TGen", "HSA BI", "HSA UPenn", "UCL BI");
seperator <-"/"
################# Input files ######################## 
# make sure to modify the paths to the correct ones
base_dir <- "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Burair_pan_scripts/Normal_tumor_pair"
meta_data_file <- paste(base_dir, "pancancer\\metadata\\merged\\metadata_summary_01_13_2021.txt", sep=seperator);
base_dir <- paste(base_dir, "pancancer\\QC\\", sep=seperator);
shared_variants_fill_list <- list("MT CUK"=paste(base_dir, "discovery\\mammary_normal_tumor_match.proportions.txt", sep=seperator),
                                  "GLM JL"=paste(base_dir, "discovery\\glioma_normal_tumor_match.proportions.txt", sep=seperator),
                                  "OM SC1"=paste(base_dir, "discovery\\melanoma_normal_tumor_match.proportions.txt", sep=seperator),
                                  "Broad PRJNA247493"=paste(base_dir, "discovery\\PRJ93_normal_tumor_match.proportions.txt", sep=seperator),
                                  "HSA BI"=paste(base_dir, "discovery\\hemangiosarcoma_normal_tumor_match.proportions.txt", sep=seperator),
                                  "MT SNU"=paste(base_dir, "validation\\PRJNA552905_MC_WES.proportions.txt", sep=seperator),
                                  #					"MT UGA"=paste(base_dir, "validation\\lab_MT_WES.proportions.txt", sep=seperator),
                                  #					"OM SC2"=paste(base_dir, "validation\\PRJEB7540_OM.proportions.txt", sep=seperator),
                                  "OSA TGen"=paste(base_dir, "validation\\Tgen_TrueN_Tpairs_WES.proportions.txt", sep=seperator),
                                  #					"OSA SC"=paste(base_dir, "validation\\PRJEB7540_OSA.proportions.txt", sep=seperator),
                                  "HSA UPenn"=paste(base_dir, "validation\\PRJNA417727_HSA_WES.proportions.txt", sep=seperator));

failed_sequencing_file <- paste(base_dir, "wgs_fail.txt", sep=seperator);

################# Output files ######################## 
# make sure to modify the paths to the correct ones
output_pdf <- paste(base_dir, "tumor_normal_pair_QC_results.pdf", sep=seperator);
output_csv <- paste(base_dir, "tumor_normal_pairing_jitterplot_data.csv", sep=seperator);
output_QC_results <- paste(base_dir, "tumor_normal_pair_QC_results.txt", sep=seperator);

################ Main code ############################

meta_data <- read.table(meta_data_file, sep="\t", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE);
meta_data <- meta_data[which(meta_data[,"Status"] == "Tumor"),];
rownames(meta_data) <- meta_data[,"SampleName"];

failed_sequencing_data <- read.table(failed_sequencing_file, row.names = 1,sep="\t", header=TRUE);
failed_sequencing_samples <- rownames(failed_sequencing_data);


# Now create the data frame for the boxplot
shared_fractions <- c();
sample_names <- c();
metric_value <- c(); # Self or Self-Best non-self
datasets <- c();
tumor_types <- c();
metric_levels <- c("Self", "Self - non-self best");

for(superset in names(shared_variants_fill_list)) {
  superset_data <- read.table(shared_variants_fill_list[[superset]], sep="\t", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1);
  failed_temp <- which(rownames(superset_data) %in% failed_sequencing_samples);
  if(length(failed_temp) > 0) {
    superset_data <- superset_data[-failed_temp, -failed_temp];
  }
  
  for(sample_index in 3:nrow(superset_data)) {
    matched_fraction <- superset_data[sample_index, sample_index]; # Either column 2 or sample_index, both are the same
    sample_name <- colnames(superset_data)[sample_index];
    cancer_type <- meta_data[sample_name, "DiseaseAcronym"];
    dataset_name <- meta_data[sample_name, "Dataset"];
    
    if(superset %in% c("OM SC1", "HSA UPenn")) {
      # These datasets contain metastasis samples, we must handle the Best non-self differently inside them
      sample_name_no_dashes <- strsplit(sample_name, split="-", fixed=TRUE)[[1]][1];
      case_sample_indices <- which(grepl(sample_name_no_dashes, colnames(superset_data), fixed=TRUE) == TRUE);
      shared_fractions_with_other_tumors <- unlist(superset_data[sample_index, -c(1, 2, case_sample_indices)]);
      shared_fractions_with_other_normals <- unlist(superset_data[-c(1, 2, case_sample_indices), sample_index]);
    } else {
      shared_fractions_with_other_tumors <- unlist(superset_data[sample_index, -c(1, 2, sample_index)]);
      shared_fractions_with_other_normals <- unlist(superset_data[-c(1, 2, sample_index), sample_index]);
    }
    
    max_other <- max(c(shared_fractions_with_other_tumors, shared_fractions_with_other_normals), na.rm=TRUE);
    if(max_other > matched_fraction) {
      best_other <- unique(c(which(superset_data[sample_index,] > matched_fraction), which(superset_data[,sample_index] > matched_fraction)));
      print(paste(dataset_name, sample_name, ", Best other:", paste(colnames(superset_data)[best_other], collapse=", ")));
    }
    shared_fractions <- c(shared_fractions, matched_fraction, (matched_fraction-max_other));
    metric_value <- c(metric_value, metric_levels);
    datasets <- c(datasets, dataset_name, dataset_name);
    tumor_types <- c(tumor_types, cancer_type, cancer_type);
    
    sample_names <- c(sample_names, sample_name);
  }
}

metric_value <- factor(metric_value, levels=metric_levels);
datasets <- factor(datasets, levels=dataset_order);
tumor_types <- factor(tumor_types, levels=cancer_types);

jitterplot_data <- data.frame("SampleName"=rep(sample_names, each=2), "DiseaseAcronym"=tumor_types, "Dataset"=datasets, "Shared_fraction"=shared_fractions, "Metric"=metric_value, check.names=F);
write.csv(jitterplot_data[order(jitterplot_data[, "Dataset"]), ], file=output_csv, row.names=FALSE, quote=FALSE);

output_QC <- data.frame("SampleName"=sample_names, "DiseaseAcronym"=meta_data[sample_names, "DiseaseAcronym"], "Dataset"=meta_data[sample_names, "Dataset"], 
                        "SelfMatch"=shared_fractions[which(metric_value == "Self")], "DiffFromBest"=shared_fractions[which(metric_value == "Self - non-self best")]);
output_QC[, "PairingQCResult"] <- "Passed";
output_QC[which(output_QC[, "DiffFromBest"] < 0), "PairingQCResult"] <- "Failed";
write.table(output_QC[order(output_QC[, "Dataset"]), ], file=output_QC_results, row.names=FALSE, sep="\t", quote=FALSE);

# Barplot data
bp_counts <- c();
bp_fractions <- c();
bp_datasets <- c();
bp_metric_valuees <- c();
bp_metric_levels <- c("Self", "Another");
bp_annt <- c();

for(bp_dataset in dataset_order) {
  sample_indices <- which(jitterplot_data[, "Dataset"] == bp_dataset);
  sample_indices <- intersect(sample_indices, which(jitterplot_data[, "Metric"] == "Self - non-self best"));
  self_count <- length(which(jitterplot_data[sample_indices, "Shared_fraction"] >= 0));
  another_count <- length(sample_indices) - self_count;
  bp_counts <- c(bp_counts, self_count, another_count);
  bp_fractions <- c(bp_fractions, self_count/length(sample_indices), another_count/length(sample_indices));
  bp_datasets <- c(bp_datasets, bp_dataset, bp_dataset);
  bp_metric_valuees <- c(bp_metric_valuees, bp_metric_levels);
  bp_annt <- c(bp_annt, paste(self_count, "/", length(sample_indices), sep=seperator), seperator);
}

bp_datasets <- factor(bp_datasets, levels=rev(dataset_order));
bp_metric_valuees <- factor(bp_metric_valuees, levels=bp_metric_levels);
barplot_data <- data.frame("counts"=bp_counts, "y"=bp_fractions, "x"=bp_datasets, "fill"=bp_metric_valuees, check.names=F);

#pdf(output_pdf, height=4.8, width=6.2);
pdf(output_pdf, height=5, width=6.84);
# plot the jitter plot
fill_colors=c("#008000", "#8B008B");
dot_size <- 2;
abs_text_size <- 20;
regular.text <- element_text(colour="black",size=abs_text_size);
xangle <- 45;
xaxis_just <- ifelse(xangle > 0, 1, 0.5);
# modify this spacing parameters between groups as desired
group_space <- 0.85;

p <- ggplot(jitterplot_data, aes(x=Dataset, y=Shared_fraction, fill=Metric)) + geom_jitter(aes(color=Metric), size=dot_size, shape=20, position=position_jitterdodge());
p <- p + facet_grid(. ~ DiseaseAcronym, scales="free_x", space="free_x");
p <- p + stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, position="dodge", geom="crossbar", size=0.2, width=0.8, color = "black", show.legend=FALSE);
p <- p + scale_fill_manual(guide=FALSE, values=fill_colors, drop=FALSE);
p <- p + scale_color_manual(name=seperator, values=fill_colors, drop=FALSE, guide=FALSE);
p <- p + labs(y="Fraction of shared variants", title="Tumor/normal pairing accuracy", fill=seperator);
p <- p + geom_hline(yintercept=0, linetype="longdash", color="yellow4", size=0.7);

p <- p + theme(axis.text=regular.text, axis.title.y=regular.text, axis.title.x=element_blank(), axis.text.x=element_text(angle=xangle, hjust=xaxis_just),
               plot.title=element_text(face="plain", size=abs_text_size, hjust=0.5),
               panel.background=element_blank(), axis.line=element_line(color="black"), strip.background=element_rect(color="black", fill="transparent", size=1.5),
               strip.text = element_text(hjust=0.5, size=abs_text_size, face="plain", color="black"), panel.spacing=unit(group_space, "lines"));

print(p);

# plot the summary
bar_plot <- ggplot(barplot_data, aes(x=barplot_data$x, y=barplot_data$y, fill=barplot_data$fill)) + 
  geom_bar(stat="identity", position=position_stack(reverse=TRUE), show.legend=TRUE, color="black", linetype="solid", alpha=0.9) + 
  coord_flip() + labs(title=seperator, x=seperator, y="Best match", fill=seperator) +
  scale_fill_brewer(palette="Set3")+
  geom_text(label=bp_annt, color="black", aes(x=x, y=y/2), size=16/(ggplot2:::.pt), hjust=0.5, fontface=1) + 
  theme_bw(base_size=abs_text_size) + theme(legend.text=element_text(size=rel(1)), legend.position="top", axis.text=element_text(size=rel(1), color="black"), axis.text.x=element_blank(),
                                            plot.background=element_blank(), axis.ticks=element_blank(), panel.border=element_blank(), 
                                            axis.line=element_blank(), panel.grid=element_blank(), plot.margin=margin(t=0.1, r=1, b=0.1, l=1, unit="in"));
print(bar_plot);
dev.off();

