# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 23:41:45 2019

@author: abc73_000
"""

with open('C:/Users/abc73_000/Desktop/java_test.txt','r') as f:
    file = f.read()
    
total = file.split('\n')[:-1]
sum = 0
for i in range(1,len(total)):
    for j in range(2,len(total[i].split('\t'))):
        each = int (total[i].split('\t')[j])
        sum+=each
    