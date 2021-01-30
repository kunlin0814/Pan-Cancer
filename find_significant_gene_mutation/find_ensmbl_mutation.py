# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:52:10 2019

@author: abc73_000
"""

"""
first argument is file path contains gene name and gene location
second argument is the gene_name file path
The output file will create a file that contains mutation number in a given gene list
ex: id1 id2 id3
    0   1   0
    

"""

import sys
import pandas as pd
candidate_info = sys.argv[1]
# r"C:\Users\abc73\Desktop\TestTMB\Pan_cancer_sign_target_ensembl_location.txt"
#
#'C:\\Users\\abc73_000\\Desktop\\Retro_gene\\new_retro_gene_list_information_CanFam3.1.99gtf.txt'
#
mut_gene_info = sys.argv[2]
#r"C:\Users\abc73\Desktop\TestTMB\CMT-2.Mutect.vcf-5steps_orientBias3-avinput.exonic_variant_function_WithGeneName"
#

sample_name = sys.argv[3]

output = open(sys.argv[4], 'w')

ensembl_mut_sum = {}
with open(candidate_info, 'r') as f:
    for each_line in f:
        if each_line[0] != "#":
            ensembl = each_line.split("\t")[0]
            ensembl_mut_sum[ensembl] = 0

mut_ensmbl_id = []

with open(mut_gene_info, 'r') as f1:
    for each_line in f1:
        ensembl = each_line.split("\t")[3].split(":")[0]
        mut_ensmbl_id.append(ensembl)


for i in mut_ensmbl_id:
    if i in ensembl_mut_sum.keys():
        ensembl_mut_sum[i] += 1

output_gene = sorted(list(ensembl_mut_sum.keys()))

# output.write("sample_names"+"\t")
# output.write("\t".join(output_gene))
# output.write("\n")
output.write(sample_name+"\t")


output_value = []
for i in output_gene:
    value = str(ensembl_mut_sum[i])
    output_value.append(value)
output.write("\t".join(output_value))
output.write("\n")


# check = pd.read_csv(r"C:\Users\abc73\Desktop\TestTMB\new_test_out.txt",sep ="\t",
#                     header =None)
