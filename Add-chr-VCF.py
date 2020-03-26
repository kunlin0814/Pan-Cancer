# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 16:46:57 2020

@author: abc73_000
"""
output = open("C:\\Users\\abc73_000\\Desktop\\final_sample.vcf", 'w')

with open("C:\\Users\\abc73_000\\Desktop\\sample.vcf", 'r')as f:
    for i in f:
        if i[0]=='#':
            pass
        else:
            info = i.split('\t')
            for j in range(len(info)):
                if j == 0:
                    output.write('chr'+str(info[j]))
                elif j==len(info)-1:
                    output.write(str(info[j]))
                else:
                    output.write(str(info[j])+'\t')

output.close()



