#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:52:10 2019

@author: abc73_000
"""

"""
first argument is retro_genefile information path contains gene name and gene location
second argument is the VCF file path

"""
import sys

retro_gene_info = sys.argv[1]
#
vcf_info = sys.argv[2]
#
with open(retro_gene_info,'r') as f:
    candidate_file = f.read()
    candidate_list =candidate_file.split('\n')[1:-1]

with open(vcf_info,'r') as f1:
    file = f1.read()
    vcf_file = file.split('\n')[:-1]

cnt = 0
for i in range(len(vcf_file)):
    chr_number = vcf_file[i].split('\t')[0]
    loc = int(vcf_file[i].split('\t')[1])
    for j  in range(len(candidate_list)):
        other_chr = candidate_list[j].split('\t')[1]
        if chr_number == other_chr :
            first_pos = int(candidate_list[j].split('\t')[2])
            last_pos = int(candidate_list[j].split('\t')[3])
            if (first_pos<=loc<= last_pos) :
                cnt+=1
       

print(cnt)




