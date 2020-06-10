# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 20:05:36 2020

@author: abc73_000
"""

import sys
import re
### create a Dict of Alt allele frequency from tumor samples
### Before five steps filtering

"""
file = []
with open("C:\\Users\\abc73_000\\Desktop\\Mutect1_Test\\Mutect2PASS.txt",'r')as f:
    for i in f :
        if (i[0] == '#'):
            pass
        else :
            file.append(i)
"""

file0 = sys.argv[1]
#'C:\\Users\\abc73_000\\Desktop\\Mutect1_Test\\Mutect2PASS.txt'

#'sys.argv[1]'

pass0 = sys.argv[2]
#'C:\\Users\\abc73_000\\Desktop\\Mutect1_Test\\Mutect2PASS.txt'



with open(file0,'r') as f:
	file = f.read()
    
    
with open(pass0,'r') as f:
	pas = f.read()
	
lst = file.split('\n')[:-1]
lst2 = pas.split('\n')[:-1]
	# output
    
out = open(pass0 + '_filteredMut','w')

dic = {}
failed = []
for i in range(len(lst)):
    normal_info = lst[i].split('\t')[10]
    tumor_info = lst[i].split('\t')[9]
    chrom = lst[i].split('\t')[0]
    pos = lst[i].split('\t')[1]
    tRef = int(tumor_info.split(':')[1].split(',')[0])
    tAlt = int(tumor_info.split(':')[1].split(',')[1])
    nRef = int(normal_info.split(':')[1].split(',')[0])
    nAlt = int(normal_info.split(':')[1].split(',')[1])
    total_tumor_depth = tRef + tAlt # total tumor read depth
    
    vaf = float(tAlt) / total_tumor_depth # VAF for tumor
    if total_tumor_depth < 10:
        failed.append(i)
    else:
        if vaf < 0.05:
            failed.append(i)
        else:
            if tAlt <= 5 and vaf < 0.15:
                failed.append(i)
            else:
                if total_tumor_depth < 20 and vaf < 0.2:
                    failed.append(i)
                else:
                    if nAlt >= 3:
                        failed.append(i)
                    else:
                        if chrom not in dic.keys():
                            dic[chrom] = []
                        dic[chrom].append(pos)
                        
                        
                        
for i in range(len(lst2)):
    info = lst2[i].split('\t')
    chrom = info[0]
    pos = info[1]
    if chrom in dic.keys():
        if pos in dic[chrom]:
            string = lst2[i] + '\n'
            out.write(string)
out.close()

    
    
#test = after5steps[2].split('\t')[9]
#re.search(r'(\d+:\d+,\d+:)(.:\d+:)(.)+',test).group(0)
 