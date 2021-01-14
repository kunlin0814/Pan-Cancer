meta_data <- read_excel("G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Burair_pan_scripts\\breed_prediction_test\\Pan-Cancer-Breed_prediction\\Breed_cluster.xlsx",
                   sheet = "Merge_InWGS")

meta_data <- setDT(meta_data)
seperator <- "/"
file_base_dir <- #"/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Burair_pan_scripts/breed_prediction_test/Pan-Cancer-Breed_prediction"
  "G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/Burair_pan_scripts/breed_prediction_test/Pan-Cancer-Breed_prediction"

Breed_QC <- function(x,y){
  if (x == "NA"){
    return("NA")
  }else if(x==y){
    return("PassBreed")
  }
  else{
    return("FailBreed")
  }
  
}

summary <- NULL
for (i in 1:nrow(meta_data)){
  x <- meta_data$Breed[i]
  y <- meta_data$BreedCluster[i]
  
  value <- Breed_QC(x,y)
  summary <- c(summary,value)
}

write.table(summary, "C:\\Users\\abc73\\Desktop\\breed.txt",col.names = F, row.names = F, sep = "\n",quote = F)

meta_data$BreedQC <- summary

final_breed <- meta_data[, "Breed"]; 
predicted_breed_indices <- intersect(which(meta_data[,.(Breed)] =="NA"), which(meta_data[, "BreedCluster"] != "Unknown"));
final_breed[predicted_breed_indices] <- meta_data[predicted_breed_indices, "BreedCluster"];
final_breed[which(meta_data[, "BreedQC"] == "FailBreed")] <- NA;
meta_data[, "FinalBreed"] <- final_breed;


write.table(meta_data, "C:\\Users\\abc73\\Desktop\\final_breed_label.txt",
            col.names = T, row.names = F, sep = "\t",quote = F,
            na = "NA")

which(meta_data[,.(Breed)] =="NA")
