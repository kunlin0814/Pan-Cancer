# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 16:04:22 2020

@author: abc73_000
"""
import sys

gene_list_file = sys.argv[1]

#'C:/Users/abc73_000/Desktop/Unique_newMC_TotalMut.txt'
MutationData= sys.argv[2]
#'C:/Users/abc73_000/Desktop/004_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut-avinput.exonic_variant_function_WithGeneName'
SampleName = sys.argv[3]
#

GeneList = []

with open (gene_list_file,'r') as f:
    file = f.read().split("\n")[:-1]
    
for gene in file:
    GeneList.append(gene)
    
GeneList.sort()        

with open (MutationData,'r') as f:
    file = f.read().split("\n")[:-1]


fileGeneList=[]
for eachLine in file :
     file_gene = eachLine.split("\t")[2]
     fileGeneList.append(file_gene)
     
MutationSummary ={}     
for gene in GeneList:
    if gene in fileGeneList:
        MutationSummary[gene]=1
    else:
        MutationSummary[gene]=0
        
output = open(SampleName+"_GeneMatrix.txt",'w')

output.write('\t')
for gene in GeneList:
    output.write(gene+'\t')

output.write("\n")

output.write(SampleName+'\t')

for gene in GeneList:
    output.write(str(MutationSummary[gene])+'\t')
    
output.write("\n")
output.close()        

