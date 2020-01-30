#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 13:38:09 2020

@author: kun-linho
"""

# for each chromosome
# min = 23840680 
# max = 123780325

import sys

## use 50 million bp as a cut off



read_length = 101
small_interval_dict ={}
large_interval_dict ={}
interval_size ={}
with open('/Volumes/Research_Data/Pan_cancer/Mapping_source/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list', 'r') as f:
    file = f.read()
    
CDS = file.split('\n')[:-1]
    
for i in range(len(CDS)):
    chrom = CDS[i].split(':')[0]
    start = int(CDS[i].split(':')[1].split('-')[0])
    end = int(CDS[i].split(':')[1].split('-')[1])
    if (start >= 50000000 and end >= 50000000):
        if chrom not in large_interval_dict.keys():
            large_interval_dict[chrom]=[(start, end)]
        else:
            large_interval_dict[chrom].append((start, end))
    else: 
         if chrom not in small_interval_dict.keys():
            small_interval_dict[chrom]=[(start, end)]
         else:
            small_interval_dict[chrom].append((start, end))
        
for i in small_interval_dict.keys():
    chrom_end = max(small_interval_dict[i])[1]
    chrom_start = min(small_interval_dict[i])[0]
    distance = chrom_end-chrom_start
    interval_size[i]= distance

#sam_file=sys.argv[1]
#output_name=sys.argv[2]
## criteria is Reads start site <=Exon end site and exon start site <= reads end site
    
unique = 0 #3
duplicate = 0 #3
Onemapped = 0 #5,9
incorrect = 0 #1
unmapped = 0 #13
with open('/Volumes/Research_Data/Pan_cancer/test.sam','r') as f:
    file = f.read()
    sam_file = file.split('\n')[40:-1]
    
reads_name = sam_file[0].split('\t')[0]    
reads_chr= sam_file[0].split('\t')[2]
reads_position = int(sam_file[0].split('\t')[3])
    
    
for line in f:
    file_lst = line.split('\t')
    if '@' in file_lst[0]:
        pass
    else :
        status = int(file_lst[1])%16
        if status == 5 or status == 9:
            Onemapped += 1
        elif status == 1:
            incorrect += 1
        elif status == 13:
            unmapped += 1
        elif status == 3:
            for ele in file_lst:
                if 'XT:' in ele:
                    status2 = ele.split(':')[2]
                    if status2 == 'U' or status2 == 'M':
                        unique += 1
                    elif status2 == 'R':
                        duplicate += 1

total = unique + duplicate + Onemapped + incorrect + unmapped
pairs = total / 2
