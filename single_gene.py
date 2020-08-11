#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 17:35:17 2019

@author: kun-linho
"""

import sys
file_name=sys.argv[1]
total_gene_list=sys.argv[2]
with open (total_gene_list,'r') as f:
    file=f.read()
#path={}
Score={}
#pathway_list=[]
gene_name=[]
total=file.split('\n')[:-1]
for i in range(len(total)):
    gene_info= total[i].split('\n')[0]
    #pathway_name = info[0].split('\t')[0]
    #pathway_list.append(pathway_name)
    Score[gene_info]=0
    #path[pathway_name]=gene_name
    

with open (file_name,'r') as f1:
     file1_info=f1.read()
total1=file1_info.split('\n')[:-1]
gene_mut=[]
for i in range(len(total1)):
    gene_mut.append(total1[i].split('\t')[2])
output=open('single_gene_mut_summary_'+file_name+'.txt','w')

for i in gene_mut:
    if i in Score.keys():
        Score[i]+=1        

#for i in Score.keys():
#    output.write(i+'\t')
#output.write('\n')

for j in Score.values():
    output.write(str(j)+'\t')

output.write('\n')

output.close()
