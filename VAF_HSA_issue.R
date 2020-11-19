#source(paste(lab_path, "pancancer\\r-code\\util.functions.R", sep=""));

residue_col_count <- 7;
residue_name_col_count <- 6;
analysis_type <- "synExcluded"; # synExcluded or synIncluded

base_dir <- paste("C:\\Users\\abc73_000\\Desktop", sep="");

#mutation_matrix_file <- paste(base_dir, "PanCancer_somaticMutationMatrix_", analysis_type, "_no_filter.txt.gz", sep="");
mutation_matrix_file <- paste(base_dir, "PanCancer_somaticMutationMatrix_synIncluded_rates_no_filter.txt.gz", sep="\\");
mutation_matrix_file_VAF <- paste(base_dir, "PanCancer_somaticMutationMatrix_synIncluded_rates_no_filter.txt.gz", sep="\\");
output_file <- paste(base_dir, "plots\\additional\\VAF_HSA_", analysis_type, ".pdf", sep="\\");

### building meta_data data frame
#source(paste(lab_path, "pancancer\\r-code\\build_germline_meta_data.R", sep=""));

# now make the mutation matrix
mutation_data <- read.table(mutation_matrix_file, header=F, sep="\t", check.names=F, stringsAsFactors=F);
sample_names <- mutation_data[1, -c(1:residue_col_count)];
gene_names <- apply(mutation_data[-c(1), 1:residue_name_col_count], MARGIN=1, function(x) {paste(x, collapse="_")});
names(gene_names) <- NULL;
mutation_matrix <- data.matrix(mutation_data[-c(1), -c(1:residue_col_count)]);
rownames(mutation_matrix) <- gene_names;
colnames(mutation_matrix) <- sample_names;
mutation_matrix <- mutation_matrix[, rownames(meta_data)];

# now make the mutation matrix
mutation_data_VAF <- read.table(mutation_matrix_file_VAF, header=F, sep="\t", check.names=F, stringsAsFactors=F);
sample_names <- mutation_data_VAF[1, -c(1:residue_col_count)];
gene_names <- apply(mutation_data_VAF[-c(1), 1:residue_name_col_count], MARGIN=1, function(x) {paste(x, collapse="_")});
names(gene_names) <- NULL;
mutation_matrix_VAF <- data.matrix(mutation_data_VAF[-c(1), -c(1:residue_col_count)]);
rownames(mutation_matrix_VAF) <- gene_names;
colnames(mutation_matrix_VAF) <- sample_names;
mutation_matrix_VAF <- mutation_matrix_VAF[, rownames(meta_data)];

outlier_hsa <- c("HSA_4", "HSA_8", "HSA_20");
random_hsa <- sample(setdiff(as.vector(unlist(sample_names[which(grepl("HSA_", sample_names) == TRUE)])), outlier_hsa), length(outlier_hsa));

convert_mutation_type <- function(x) {
	if(x == "A>C") {
		return("T>G");
	} else if(x == "A>G") {
		return("T>C");
	} else if(x == "A>T") {
		return("T>A");
	} else if(x == "G>A") {
		return("C>T");
	} else if(x == "G>C") {
		return("C>G");
	} else if(x == "G>T") {
		return("C>A");
	} else {
		return(x);
	}
}

count_mutation_types <- function(x, signature_levels) {
	result <- rep(0, length(signature_levels));
	names(result) <- signature_levels;
	get_table <- table(x);
	for(i in 1:length(get_table)) {
		result[names(get_table)[i]] <- get_table[i];
	}
	return(result/length(x));
}
signature_colors <- c("cyan","black","red","gray","green","pink");
signature_levels <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G");

analyzed_samples <- c(outlier_hsa, random_hsa);
y <- c();
x <- c();
signatures_y <- c();
signatures_x <- c();
for(sample in analyzed_samples) {
	sample_mutations <- intersect(rownames(mutation_matrix_VAF), rownames(mutation_matrix)[which(mutation_matrix[, sample] == 1)]);
	y <- c(y, mutation_matrix_VAF[sample_mutations, sample]);
	x <- c(x, rep(sample, length(sample_mutations)));
	mutation_types <- sapply(sample_mutations, function(x) {convert_mutation_type(paste(strsplit(x, split="_", fixed=T)[[1]][4:5], collapse=">"))});
	signatures_y <- c(signatures_y, count_mutation_types(mutation_types, signature_levels));
	signatures_x <- c(signatures_x, rep(sample, length(signature_levels)));
}

count_mutation_types()

x <- factor(x, analyzed_samples);
signatures_x <- factor(signatures_x, analyzed_samples);
signatures_fill <- factor(names(signatures_y), signature_levels);

xangle <- 0;

my_jitter <- function(dataframe, dot_size=1, abs_text_size=20, xangle=30) {
	regular.text <- element_text(colour="black",size=abs_text_size);
	xaxis_just <- ifelse(xangle > 0, 1, 0.5);
	
	p <- ggplot(dataframe, aes(x=x, y=y)) + geom_jitter(color="black", size=dot_size, shape=20, position=position_jitter(width=0.2));
	p <- p + stat_summary(fun.y=median, fun.ymin=median, fun.ymax=median, geom="crossbar", size=0.2, width=0.8, color = "red", show.legend=FALSE);
	p <- p + theme(axis.text=regular.text, axis.title=regular.text, axis.text.x = element_text(angle=xangle, hjust=xaxis_just), legend.position="none",
				panel.background=element_blank(), axis.line=element_line(color="black"), 
				plot.margin=margin(1,1,1,3, unit="cm"));
	return(p);
}

pdf(output_file, height=6.5, width=10);
	data <- data.frame(x=x, y=y);
	p <- my_jitter(data, xangle=xangle);
	p <- p + labs(title="", y="VAF", x="");
	plot(p);
	
	data <- data.frame(x=signatures_x, y=signatures_y, fill=signatures_fill);
	p <- ggplot(data, aes(x=x, y=y, fill=fill)) + geom_bar(stat="identity", position="stack") + scale_fill_manual(values=signature_colors);
	plot(p);
dev.off();

