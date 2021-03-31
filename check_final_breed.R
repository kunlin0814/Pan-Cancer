library(data.table)

original <- fread("/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/final_breed_data_02_18.txt")
new <- fread("/Volumes/Research/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/arrange_table/Final_breed_labels_03_30.txt")
new_breed <- new[DiseaseAcronym=="OM",.(SampleName,FinalBreed)]
new_om <- new_breed[grepl("DD",new_breed$SampleName),]
original_breed <- original[grepl("DD",original$Case_ID),]
final_result <- merge(new_om,original_breed,by.x ='SampleName',by.y = 'Case_ID'
                      ,all.y = T)
original_breed$Case_ID


