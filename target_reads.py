#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 13:38:09 2020

@author: kun-linho
"""

# for each chromosome
# min = 23840680 
# max = 123780325

## use 50 million bp as a cut off
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
