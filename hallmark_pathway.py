#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 14:13:04 2019

@author: kun-linho
"""
import sys
file_name=sys.argv[1]
total_pathway_genes=sys.argv[2] # hallmark gene names
with open (total_pathway_genes,'r') as f:
    file=f.read()
total=file.split('\n')[:-1]
path={}
pathway_list=[]
score={}
for i in range(len(total)):
    info= total[i].split('\n')
    pathway_name = info[0].split('\t')[0]
    pathway_list.append(pathway_name)
    gene_name=info[0].split('\t')[3:-1]
    path[pathway_name]=gene_name
    score[pathway_name]=0
    #globals(pathway_name)
    
with open (file_name,'r') as f1:
    file1_info=f1.read()
total1=file1_info.split('\n')[:-1]
gene_mut=[]
for i in range(len(total1)):
    gene_mut.append(total1[i].split('\t')[2])
output=open('hallmark_gene_pathway_mut_summary_'+file_name+'.txt','w')

for i in gene_mut:
    for j in score.keys():
        if i in path[j]:
            score[j]+=1
        
for i in score.keys():
    output.write(i+'\t')
output.write('\n')

for j in score.values():
    output.write(str(j)+'\t')

output.write('\n')

output.close()