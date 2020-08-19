# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 23:19:14 2020

@author: abc73_000
"""
import sys

annovarFile= sys.argv[1]
#"C:/Users/abc73_000/Desktop/GermlineSomaticOverlap/CMT-2_GermlineSomaticSNVlOverlap.vcf-avinput.exonic_variant_function"
OverlapFile= sys.argv[2]
#"C:/Users/abc73_000/Desktop/GermlineSomaticOverlap/CMT-2_GermlineSomaticSNVlOverlap.vcf"
RetroGeneList="/work/szlab/kh31516_Lab_Share_script/Retro_gene_source/new_retro_gene_list_CanFam3.1.99gtf.txt"
#"C:/Users/abc73_000/Desktop/GermlineSomaticOverlap/new_retro_gene_list_CanFam3.1.99gtf.txt"
output=open(sys.argv[3],'w')

with open(RetroGeneList,'r') as f:
    RetroGeneList = f.read().split("\n")[:-1]

with open(annovarFile,'r') as f:
    annovar = f.read().split("\n")[:-1]

with open(OverlapFile, 'r') as f:
    Overlap = f.read().split("\n")[:-1]

for i in range(len(annovar)):     
    ensemblID= annovar[i].split("\t")[2].split(":")[0]
    if ensemblID not in RetroGeneList:
        output.write(Overlap[i])
        output.write("\n")
        

output.close()