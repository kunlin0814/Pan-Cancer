#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:06:12 2020

@author: abc73_000
"""
import os
import sys
import numpy as np

with open("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/PanCancerCoverage/MC_Coverage/CMT-2/SRR7780922_Normal_Coverage_distribution.txt", 'r')as f:
    Normal = f.read().split("\n")[:-1]

with open("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/PanCancerCoverage/MC_Coverage/CMT-2/SRR7780923_Tumor_Coverage_distribution.txt", 'r')as f:
    Tumor = f.read().split("\n")[:-1]

output = open("/Volumes/Research_Data/MAC_Research_Data/Pan_cancer/Pan_cancer-analysis/PanCancerCoverage/MC_Coverage/CMT-2/testlog2.txt",'w')

Normallocation = []
normalCoverage=[]
    
for eachLine  in Normal:
    pos = eachLine.split("\t")[0]
    coverage = float(eachLine.split("\t")[1])
    normalCoverage.append(coverage)
    Normallocation.append(pos)

normalCoverage = np.array(normalCoverage)
normalMean = np.mean(normalCoverage)
normalStd = np.std(normalCoverage)

normalizedNormal = (normalCoverage-np.min(normalCoverage))/(np.max(normalCoverage)-np.min(normalCoverage))
#normalCoverage/normalMean
#(normalCoverage-np.min(normalCoverage))/(np.max(normalCoverage)-np.min(normalCoverage))
    
        


tumorCoverage = []

for eachLine in Tumor:
    pos = eachLine.split("\t")[0]
    coverage = float(eachLine.split("\t")[1])
    tumorCoverage.append(coverage)
    #tumorLocation.append(pos)
    


tumorCoverage = np.array(tumorCoverage)
tumorMean = np.mean(tumorCoverage)
tumorStd = np.std(tumorCoverage)

normalizedTumor = (tumorCoverage-np.min(tumorCoverage))/(np.max(tumorCoverage)-np.min(tumorCoverage))
#normalizedTumor/tumorMean
#(tumorCoverage-np.min(tumorCoverage))/(np.max(tumorCoverage)-np.min(tumorCoverage))


normalizedNormal = normalizedNormal+0.001
normalizedTumor = normalizedTumor+0.001


log2Ration = np.log2(normalizedNormal/normalizedTumor)

for i in range(len(Normallocation)):
    interval = Normallocation[i]
    log2 = str(log2Ration[i])
    output.write(interval+"\t"+log2+"\n")
    
output.close()








