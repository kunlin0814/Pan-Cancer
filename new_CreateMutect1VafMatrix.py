#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Mon Jul 27 19:30:21 2020

@author: abc73_000
"""
### the script can get the mutect VAF info

### the script need 3 input files
# 1: total gene mutation and it's location, 
# 2: Annovar result with gene name (from Mutect1, then annovar)
# 3: Mutect Stats result (to get the alleles info)

## other two parameters 
# 4 : outputfile_name
# 5 : sample_name

import sys

file0 = sys.argv[1]
#
#'C:\\Users\\abc73_000\\Desktop\\uniq_total_ProteinMutationResults.txt'
#
sampleNonsyn = sys.argv[2]
#
#'C:\\Users\\abc73_000\\Desktop\\NonSyn_CMT-2_Before_filtering_WithGeneName'

vafFile = sys.argv[3]
#
#"C:\\Users\\abc73_000\\Desktop\\CMT-2_PASS.stat"

sampleName = sys.argv[4]
#'CMT-2'
#


### Create Mutation Dict from Total Mutation Database ###
with open(file0,'r') as f:
    file = f.read().split("\n")[:-1]
    
Mutation_Data = {}

for each_line in file :
    chrom = each_line.split("\t")[0]
    start = each_line.split("\t")[1]
    #end = each_line.split("\t")[2]
    mutation = each_line.split("\t")[3]
    key = chrom+":"+start
    if key not in Mutation_Data.keys():
        Mutation_Data[key] = [mutation]
    else:
        Mutation_Data[key].append(mutation)

Mutation_Data_keys= list(Mutation_Data.keys())
Mutation_Data_keys.sort()

##### Create a Mutation summary from the sample

Sample_Dict = {}
with open(sampleNonsyn, 'r')as f:
    sample = f.read().split('\n')[:-1]

for each_line in sample:
    gene_name = each_line.split("\t")[2]
    chrom = each_line.split("\t")[4]
    start= each_line.split("\t")[5]
    #end = each_line.split("\t")[6]
    mutation = each_line.split("\t")[3].split(":")[-1][2:-1]
    value = gene_name+'_'+mutation
    key = chrom+":"+start
    if key not in Sample_Dict.keys():
        Sample_Dict[key] = value

        
### Get the VAF info created from Mutect Stat    
    
with open(vafFile,'r') as f:
    vaf_info = f.read().split('\n')[:-1]

VAF_summary={}        
for each_line in vaf_info:
    info = each_line.split('\t')
    chrom = info[0]
    pos = info[1]
    tRef = int(info[2])
    tAlt = int(info[3])
    nRef = int(info[4])
    nAlt = int(info[5])
    total_tumor_depth = tRef + tAlt # total tumor read depth
    vaf = float(tAlt) / total_tumor_depth # VAF for tumor
    key = chrom+":"+pos
    value = vaf
    VAF_summary[key]= value
    
output= open(sampleName+'_Mutation_VAF_info.txt','w')

Total_summary ={}
for position in Mutation_Data_keys:
    #location= Mutation_Data[position]
    if position in Sample_Dict.keys():
        allGeneList = Mutation_Data[position]
        sampleGene = Sample_Dict[position]
        for i in allGeneList:
            if i == sampleGene:
                vafValue = VAF_summary[position]
                Total_summary[i] = vafValue
            else:
                Total_summary[i] = 0
        
    else:
        gene_list = Mutation_Data[position]
        for i in gene_list:
            Total_summary[i] = 0

"""        
output.write("\t")

for position in Mutation_Data_keys:
    geneList = Mutation_Data[position]
    for gene in geneList:
        output.write(gene+'\t')

output.write("\n")
"""
output.write(sampleName+"\t")

for position in Mutation_Data_keys:
    geneList = Mutation_Data[position]
    for gene in geneList:
        output.write(str(Total_summary[gene])+'\t')

output.write("\n")    
        
output.close()    
    


