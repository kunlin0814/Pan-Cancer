#!/usr/bin/python

import sys
file_name=sys.argv[1]
total_pathway_genes=sys.argv[2]

#file_name='/Users/kun-linho/Desktop/Test_DNA_repair/non_syn_CMT-262_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS-avinput.exonic_variant_function_WithGeneName'
#total_pathway_genes='/Users/kun-linho/Desktop/Test_DNA_repair/DNA_repair_gene.txt'
with open (total_pathway_genes,'r') as f:
    file=f.read()
path={}
pathway_list=[]
score={}
total=file.split('\n')[:-1]
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

output=open('DNA_repair_mut_summary_'+file_name+'.txt','w')

for i in gene_mut:
    if i in score.keys():
        score[i]=1
        #print(i)

        

for i in score.keys():
    output.write(i+'\t')
output.write('\n')

for j in score.values():
    output.write(str(j)+'\t')

output.write('\n')

output.close()
