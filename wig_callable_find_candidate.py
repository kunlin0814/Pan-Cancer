# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 17:31:19 2019
@author: abc73_000
"""
"""
The script takes two input files. The first one is wig callable base file 
, and the second takes the gene name list you want to identify the callable bases for these genes (ensemble id)

The output files will give you the callable base for each candidate genes (ensemble ID)

The function is the same as the java script "GetCallableCounts.java" but in python version

"""
with open('/Users/kun-linho/Desktop/Retro_gene_finding/CMT-2_rg_added_sorted_dedupped_removed.bam_coverage.wig.txt', 'r')as f:
    file = f.read()
    total_line = file.split('\n')[1:-1]

GeneChromDict={}
GeneLocDict ={}
Total_list = []
with open('/Volumes/Research_Data/Pan_cancer/Retro_gene_finding/RetroGeneList/retro_gene_list_information.txt','r') as f:
    candidate_file = f.read()
    candidate_list =candidate_file.split('\n')[1:-1]

for i in range(len(candidate_list)):
    id = candidate_list[i].split('\t')[0]
    chro = candidate_list[i].split('\t')[1]
    start = candidate_list[i].split('\t')[2]
    end = candidate_list[i].split('\t')[3]
    GeneChromDict[id] = chro
    GeneLocDict[id] = (start,end )
    Total_list.append(id)

CurrentIDCnt = 0
RetroTotalCnt = 0
currentStep = 0
currentpos = 0  
startedGenesInsideCurrentChr=[]
unstartedGenesInsideCurrentChr =[]
unstartedGenes = Total_list
temp = []
currentChr = ' '
result ={}

  
for i in range(len(total_line)):
    if "step" in total_line[i]:
        chrosome = total_line[i].split(' ')[1].split('chrom=')[1]
        currentpos = int(total_line[i].split(' ')[2].split('start=')[1])
        #currentStep = int(total_line[i].split(' ')[3][5:])
        if (chrosome!=currentChr):
            startedGenesInsideCurrentChr.clear()
            unstartedGenesInsideCurrentChr.clear()
            for unstartgene in unstartedGenes :
                GeneChrom =  GeneChromDict[unstartgene]
                if GeneChrom == chrosome :
                    unstartedGenesInsideCurrentChr.append(unstartgene)
            unstartedGenes= list(filter(lambda i: i not in unstartedGenesInsideCurrentChr, unstartedGenes))
            currentChr = chrosome
            print('current '+ currentChr)
    elif ((total_line[i] == '1') or (total_line[i]=='0') ):
        #print ('Start Counting')
        currentpos += 1
        for k in unstartedGenesInsideCurrentChr :
            if (int(GeneLocDict[k][0]) > int(currentpos)):
                pass
            elif ((int(GeneLocDict[k][0])<= int(currentpos)) and (int(GeneLocDict[k][1])>= currentpos)) :
                startedGenesInsideCurrentChr.append(k)
                temp.append(k)
                #print (k + ' is in the interval')
            else: 
                temp.append(k)
        unstartedGenesInsideCurrentChr = list(filter(lambda i: i not in temp, unstartedGenesInsideCurrentChr))
        temp.clear();
        for gene in startedGenesInsideCurrentChr :
            location = GeneLocDict[gene][1]
            if (int(location) < currentpos):
                temp.append(gene) 
        startedGenesInsideCurrentChr = list(filter(lambda i: i not in temp, startedGenesInsideCurrentChr))
        #Total_list = list((filter(lambda i: i not in temp, Total_list)))
        temp.clear()
        if (total_line[i]=='1'):
            RetroTotalCnt+=1
            for gene in startedGenesInsideCurrentChr:
                if gene in result.keys():
                    result[gene]+=1
                else:
                    result[gene]= 1
        #
                
    
for i in Total_list:
    if i not in result.keys():
        result[i]=0
              
output = open ('/Users/kun-linho/Desktop/finalized.txt','w')
for i in result.keys():      
    output.write(str(i)+'\t')

output.write('\n')

for i in result.keys():      
    output.write(str(result[i])+'\t')

output.write('\n')
    

output.close()

"""
sum1=0
for i in result.values():
    sum1+=i

print(sum1)
"""
        
  
"""       
        
        Total_list = filter(lambda i: i not in temp, Total_list)
        temp.clear()
    print('the len of the gene ' + str(len(list(Total_list))))
    #print(RetroTotalCnt)

for i in total_line:
    print(i)




   
start = []
for i in range(len(total_line)):
    start_pos = int(total_line[i].split(" ")[2].split('start=')[1])
    chrom = total_line[i].split(" ")[1].split('chrom=')[1]
    for j in range(len(candidate_list)) :
        if (chrom == candidate_list[j].split('\t')[1]):
            if (int(candidate_list[j].split('\t')[2])<=start_pos<=int(candidate_list[j].split('\t')[2])):
                print('wig file is')
                print(chrom)
                print(start_pos)
                print('candidate_list  is')
                print(candidate_list[j])
"""