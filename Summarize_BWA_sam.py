#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 10:53:35 2020

@author: kun-linho
"""

import sys

#sam_file=sys.argv[1]
sam_file= '/Volumes/Research_Data/Pan_cancer/Mapping_source/Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list'
#output_name=sys.argv[2]


#with open(sam_file,'r') as f:
#	file=f.read()

#summary=open(output_name+'_Mapping_summary.txt','w')

unique = 0 #3
duplicate = 0 #3
Onemapped = 0 #5,9
incorrect = 0 #1
unmapped = 0 #13
with open(sam_file,'r') as f:
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

unique_rate = float(unique)/total
dup_rate = float(duplicate)/total
Onemap_rate = float(Onemapped)/total
incorrect_rate = float(incorrect)/total
unmapped_rate = float(unmapped)/total

#summary.write(output_name + '\n')
#summary.write('ID\tTotal_pairs\tUniquely_mapped_rate\tRepeatedly_mapped_rate\t1read_mapped_rate\tIncorrectly_mapped_rate\tUnmapped_rate\tUniquely_mapped\tRepeatedly_mapped\t1read_mapped\tIncorrectly_mapped\tUnmapped\tTotal_reads\n')
#summary.write(output_name+'\t'+str(unique) + '\t' + str(duplicate) + '\t' + str(Onemapped) + '\t' + str(incorrect) + '\t' + str(unmapped) + '\t' + str(total) + '\t' + str(pairs) + '\n')
#summary.write(output_name+'\t'+str(pairs)+ '\t'+ str(unique_rate) + '\t' + str(dup_rate) + '\t' + str(Onemap_rate) + '\t' + str(incorrect_rate) + '\t' + str(unmapped_rate) +'\t'+str(unique)+'\t' + str(duplicate) + '\t' + str(Onemapped) + '\t' + str(incorrect) + '\t' + str(unmapped) + '\t' + str(total) + '\t\t\n')


#summary.close()

