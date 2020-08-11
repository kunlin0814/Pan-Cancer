#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:06:12 2020

@author: abc73_000
"""

import re
import sys
import string

input_file = sys.argv[1]
#"/Users/kun-linho/Desktop/CMT-102-DbSNPFiltering_Consensus.sINDEL.vcf"
#sys.argv[1]
output_file= sys.argv[2]
#"/Users/kun-linho/Desktop/Test.txt"
#sys.argv[2]

output = open(output_file,'w')

with open(input_file,'r')as f:
    file = f.read().split("\n")[:-1]
    
for eachLine in file:
   if eachLine.startswith('#'):
       pass
       #output.write(eachLine+"\n")
   else:
        each_line_rec = eachLine.split('\t')
        filter_status = eachLine.split('\t')[6]
        if (filter_status == "PASS" or filter_status =="LowQual"):
            #re.search(r'MVL=(\d,\d,\d)',eachLine).group(1)
            #filter_status = reviseStatus(VCFstatus)
            #each_line_rec[6] = filter_status
            #final_each_line = "\t".join(each_line_rec)    
            output.write(eachLine+"\n")

output.close()
    

