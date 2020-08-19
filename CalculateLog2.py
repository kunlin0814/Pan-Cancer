#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:06:12 2020

@author: abc73_000
"""
import os
import sys
import numpy as np

NormalFile=sys.argv[1]
TumorFile=sys.argv[2]
#"G:/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/CoverageTest/MC_Coverage/CMT-2/SRR7780922_Normal_Coverage_distribution.txt"
with open(NormalFile, 'r')as f:
    Normal = f.read().split("\n")[:-1]

with open(TumorFile, 'r')as f:
    Tumor = f.read().split("\n")[:-1]

output = open(sys.argv[3],'w')

Normallocation = []
normalCoverage=[]
    
for eachLine  in Normal:
    pos = eachLine.split("\t")[0]
    coverage = float(eachLine.split("\t")[1])
    normalCoverage.append(coverage)
    Normallocation.append(pos)

normalCoverage = np.array(normalCoverage)
normalMean = np.mean(normalCoverage)

normalizedNormal = normalCoverage/normalMean
#
#(normalCoverage-np.min(normalCoverage))/(np.max(normalCoverage)-np.min(normalCoverage))
    
        


tumorCoverage = []

for eachLine in Tumor:
    pos = eachLine.split("\t")[0]
    coverage = float(eachLine.split("\t")[1])
    tumorCoverage.append(coverage)
    #tumorLocation.append(pos)
    


tumorCoverage = np.array(tumorCoverage)
tumorMean = np.mean(tumorCoverage)


normalizedTumor = tumorCoverage/tumorMean
#
#(tumorCoverage-np.min(tumorCoverage))/(np.max(tumorCoverage)-np.min(tumorCoverage))


normalizedNormal = normalizedNormal+0.001
normalizedTumor = normalizedTumor+0.001


log2Ration = np.log2(normalizedTumor/normalizedNormal)

for i in range(len(Normallocation)):
    interval = Normallocation[i]
    log2 = str(log2Ration[i])
    output.write(interval+"\t"+log2+"\n")
    
output.close()








