# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 19:30:21 2020

@author: abc73_000
"""

import sys
annovar_file = sys.argv[1]
#"C:\\Users\\abc73_000\\Desktop\\NonSyn_CMT-2_Before_filtering_WithGeneName"
#
output=open(sys.argv[2],'w')
#"C:\\Users\\abc73_000\\Desktop\\result.txt"
#sys.argv[2]

## protein Mutation info


with open(annovar_file,'r') as f:
    file = f.readlines()

for each_line in file:
    gene_name = each_line.split("\t")[2]
    chrom = each_line.split("\t")[4]
    start= each_line.split("\t")[5]
    end = each_line.split("\t")[6]
    mutation = each_line.split("\t")[3].split(":")[-1][2:-1]
    output.write(chrom+'\t'+start+'\t'+end+'\t'+gene_name+'\n')
    
output.close()