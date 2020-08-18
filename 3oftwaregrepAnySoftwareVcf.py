#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 19 17:06:12 2020

@author: abc73_000
"""
import os
import re
import sys
import string

Consenus = sys.argv[1]
#"/Users/kun-linho/Desktop/CMT-102-DbSNPFiltering_Consensus.sINDEL.vcf"
Germline = sys.argv[2]
#sys.argv[1]
output = open(sys.argv[3],'w')
#"/Users/kun-linho/Desktop/Test.txt"
#sys.argv[2]
if os.stat(Consenus).st_size == 0:
    output.write("Consensus Mutation file is empty")
else:
    with open(Consenus, 'r') as f:
        ConsensusMutation= f.read().split("\n")[:-1]
    ConsensusDict={}

    for eachLine in ConsensusMutation:
        if eachLine.startswith('#'):
            pass
        else:
            chrom = eachLine.split("\t")[0]
            pos = eachLine.split("\t")[1]
            ref = eachLine.split("\t")[3]
            alt = eachLine.split("\t")[4]
            ConsensusKey = chrom+":"+ pos
            if ConsensusKey not in ConsensusDict.keys():
                ConsensusDict[ConsensusKey] = [(ref,alt)]
            else:
                ConsensusDict[ConsensusKey].append[(ref,alt)]
            
            
            
    with open(Germline,'r') as f:
        GermlineMutation = f.read().split("\n")[:-1]
        
    GermlineDict = {}   
    for eachLine in GermlineMutation:
        chrom = eachLine.split("\t")[0]
        pos = eachLine.split("\t")[1]
        ref = eachLine.split("\t")[3]
        alt = eachLine.split("\t")[4]
        GermlineKey = chrom+":"+ pos
        if GermlineKey in ConsensusDict.keys():
            if (ref,alt) in ConsensusDict[GermlineKey]:
                output.write(eachLine)
                output.write("\n")
            
            
output.close()    
    

    

