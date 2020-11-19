# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 20:05:36 2020

@author: abc73_000
"""
import sys
import re
### create a Dict of Alt allele frequency from tumor samples
### Before five steps filtering

before5_VAF_info = "C:\\Users\\abc73_000\\Desktop\\Mutect1_test\\HSA_3_PASS.stat"
#sys.argv[1]
after5stepsVCF= "C:\\Users\\abc73_000\\Desktop\\Mutect1_test\\HSA_3_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS_filteredMut"
#sys.argv[2]
after5stepsVAF_info = "C:\\Users\\abc73_000\\Desktop\\Mutect1_test\\test_after.txt"
#sys.argv[3]

after_filtering = open(after5stepsVAF_info,'w')

allelDict={}

with open(before5_VAF_info,'r')as f:
    file = f.readlines()

for i in file :
    chrom = i.split('\t')[0]
    pos = i.split('\t')[1]
    freq = i.split('\t')[4]
    allelDictKey = str(chrom+pos)
    if allelDictKey not in allelDict.keys():
        allelDict[allelDictKey] = freq
    
#### compare with 5 steps filtering ####

with open(after5stepsVCF,'r') as f:
    after5steps = f.readlines()
    
    
for i in after5steps:
    chrom = i.split('\t')[0]
    pos = i.split('\t')[1]
    vcf_key = str(chrom+pos)
    #ad_info
    if vcf_key in allelDict.keys():
        #print(allelDict[vcf_key])
        after_filtering.write(chrom+'\t'+ str(pos)+'\t')
        after_filtering.write(allelDict[vcf_key])
        after_filtering.write("\n")
    
        
after_filtering.close()    
