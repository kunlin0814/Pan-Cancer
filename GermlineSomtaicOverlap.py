# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 19:30:43 2020

@author: abc73_000
"""


import sys

Consenus = sys.argv[1]
#'C:/Users/abc73_000/Desktop/Test_Area/SomaticSeq/AtLeastTwo_CMT-100-DbSNPFiltering_Consensus.sSNV.vcf'



Germline = sys.argv[2]
#'C:/Users/abc73_000/Desktop/Test_Area/SRR7779579-rg_added_sorted_dedupped_split.realigned.bam.filter.vcf-PASS'
#sys.argv[2]

output = open(sys.argv[3],'w')

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NORMAL	TUMOR

with open(Consenus, 'r') as f:
    ConsensusMutation= f.read().split("\n")[:-1]

ConsensusDict={}

for eachLine in ConsensusMutation:
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
        if (ref, alf) in ConsensusDict[GermlineKey]:
            output.write(eachLine)
            
            
output.close()    
    
    


