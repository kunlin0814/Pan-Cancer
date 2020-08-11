#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 14:15:45 2019

@author: kun-linho
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


target_chr = ['chr'+str(i) for i in range(1,39)]
target_chr.append('chrx')

gene_location={}
gene_transcript={}
transcript_number={}
real_gene={}
candidate_ensembleID_number = []
retro_candidate_gene=[]
total_gene={}
final_retro_gene = []
final_retro_gene_info ={}

distance_cutoff = 0
with open("G:\\Pan_cancer\\Mapping_source\\Canis_familiaris.CanFam3.1.99.chr.gtf") as f:
    file = f.read()
    
total = file.split('\n')[5:-1]

for i in range(len(total)):
    if('gene_name' in total[i] and total[i].split('\t')[2] == "gene" and 'protein_coding' in total[i]):
        ensemble_id = total[i].split('gene_id')[1].split(';')[0].replace('"', '' ).replace(" ","")
        real_gene[ensemble_id]=total[i].split('gene_id')[1].split(';')[2].split(" ")[2].replace('"', '' ).replace(" ","")
        
    elif (total[i].split('\t')[2] == "transcript"and 'protein_coding' in total[i]):
        ensemble_id = total[i].split('gene_id')[1].split(';')[0].replace('"', '' ).replace(" ","")
        if ensemble_id in gene_transcript.keys():
            gene_transcript[ensemble_id].extend([total[i].split('gene_id')[1].split(';')[2].split(" ")[2].replace('"', '' ).replace(" ","")])
        else :
            gene_transcript[ensemble_id] = [total[i].split('gene_id')[1].split(';')[2].split(" ")[2].replace('"', '' ).replace(" ","")]
    elif  total[i].split('\t')[2] == "exon":
        transcript_id= total[i].split('gene_id')[1].split(';')[2].split(" ")[2].replace('"', '' ).replace(" ","")
        if transcript_id in transcript_number.keys():
            transcript_number[transcript_id]+=1
        else :
            transcript_number[transcript_id]= 1
    elif ('gene_name' not in total[i] and total[i].split('\t')[2] == "gene" and 'protein_coding' in total[i]):
        ensemble_id = total[i].split('gene_id')[1].split(';')[0].replace('"', '' ).replace(" ","")
        retro_candidate_gene.append(ensemble_id)
        gene_location[ensemble_id] = (int(total[i].split('\t')[3]),int(total[i].split('\t')[4]))
   
     
for i in range(len(total)):
    if (total[i].split('\t')[2] == "gene"):    
        ensemble_id = total[i].split('gene_id')[1].split(';')[0].replace('"', '' ).replace(" ","")
        chromosome= total[i].split('\t')[0].split(';')[0].replace('"', '' ).replace(" ","")
        total_gene[ensemble_id]= chromosome



for i in retro_candidate_gene:
    transcript = gene_transcript[i]
    if len(transcript)==1:
        exon_number = int(transcript_number[transcript[0]])
        if exon_number==1:
            candidate_ensembleID_number.append(i)



    
gene_distance={}
for i in candidate_ensembleID_number:
    gene_distance[i]=int(gene_location[i][1]-gene_location[i][0])
    if gene_distance[i] > distance_cutoff :
       chro =  total_gene[i]
       locatio = gene_location[i]
       final_retro_gene_info[i]=[chro,locatio]


non_retro = {}          
for i in total_gene.keys():
    if i not in final_retro_gene_info.keys():
        non_retro[i] = total_gene[i]

non_retro_list = open('G:\\Pan_cancer\\Mapping_source\\new_non_retro_gene_list.txt','w')
#non_retro_list.write('Ensemble_id'+'\t'+'chromosome_location'+'\n')
for i in real_gene.keys():  
    if total_gene[i] in target_chr:
        non_retro_list.write(i+'\n')
        #non_retro_list.write(total_gene[i]+'\n')
    
non_retro_list.close()    
   
"""
plt.figure(figsize=(16,9))
sns.set(font_scale=3)    
sns.distplot(list(gene_distance.values()),kde=False, axlabel= 'Gene_length', color='black') 
#plt.xlim(0, 2000)
#plt.ylim(0,125)
plt.savefig('/Volumes/Research_Data/Pan_cancer/Retro_gene_finding/candidate_retro_gene_length.png')
plt.close()            
"""
    

output_genelist = open('G:\\Pan_cancer\\Mapping_source\\new_retro_gene_list.txt','w')
output = open('G:\\Pan_cancer\\Mapping_source\\new_retro_gene_list_information.txt','w')    
output.write('Ensemble_id'+'\t'+'chromosome_location'+'\t'+'gene_start'+'\t'+'gene_end'+'\n')
for i in final_retro_gene_info.keys():
    output_genelist.write(i+'\n')
    output.write(i+'\t')
    output.write(str(final_retro_gene_info[i][0]+'\t'))
    output.write(str(final_retro_gene_info[i][1][0])+'\t')
    output.write(str(final_retro_gene_info[i][1][1])+'\n')

output.close()
output_genelist.close()

"""    
retro_data= pd.read_csv("/Volumes/Research_Data/Pan_cancer/Retro_gene_finding/retro_gene_list_information.txt",sep='\t')
retro_data['Difference']= retro_data['gene_end']-retro_data['gene_start']
#retro_data.sort_values(by = 'chromosome_location',inplace= True)
retro_data[retro_data['Difference']>=800]
"""    
        