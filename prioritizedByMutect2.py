# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:06:12 2020

@author: abc73_000
"""

import re
import sys

def reviseStatus(VCFstatus):
    filter_status = ''
    if VCFstatus == '1,1,1':
        filter_status = 'PASS-1'
    elif VCFstatus == '1,0,1':
        filter_status = 'PASS-2-1'
    elif VCFstatus == '1,1,0':
        filter_status = 'PASS-2-2'
    elif VCFstatus == '1,0,0':
        filter_status = 'PASS-2'
    elif VCFstatus == '0,1,1':
        filter_status = 'PASS-3'
    else:
        filter_status = 'LowQual'
    return filter_status

input_file = sys.argv[1]
output_file= sys.argv[2]

output = open(output_file,'w')

with open(input_file,'r')as f:
    file = f.readlines()
    
for i in range(len(file)):
   if file[i][0] =='#':
       #output.write(file[i])
       pass
   else:
        each_line_rec = file[i].split('\t')
        filter_status = file[i].split('\t')[6]
        VCFstatus = re.search(r'MVL=(\d,\d,\d)',file[i]).group(1)
        filter_status = reviseStatus(VCFstatus)
        each_line_rec[6] = filter_status
        final_each_line = "\t".join(each_line_rec)    
        output.write(final_each_line)

output.close()
    
    

