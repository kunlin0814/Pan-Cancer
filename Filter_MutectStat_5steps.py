#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

## The script will use Mutect1 PASS file and Stat file create by grep "keep" Stat file from Mutect1 as input
## output, five steps filter Mutect1 and VAF before and AFter and also the reasons that filtering out


file0 = sys.argv[1]
#"C:\\Users\\abc73_000\\Desktop\\Mutect1_test\\HSA_3_PASS.stat"
#
pass0 = sys.argv[2]
#"C:\\Users\\abc73_000\\Desktop\\Mutect1_test\\HSA_3_rg_added_sorted_dedupped_removed.MuTect.vcf-PASS"
vaf_before = sys.argv[3]
#"C:\\Users\\abc73_000\\Desktop\\Mutect1_test\\Test_before.txt"
vaf_after = sys.argv[4]
#"C:\\Users\\abc73_000\\Desktop\\Mutect1_test\\Test_ar.txt"
#sys.argv[2]
whyout = sys.argv[5]
# "C:\\Users\\abc73_000\\Desktop\\Mutect1_test\\Test_why.txt"
#sys.argv[3]
vafbeforeout = open(vaf_before,'w')
vafafterout = open(vaf_after,'w')
with open(file0,'r') as f:
    file = f.read()
with open(pass0,'r') as f:
    pas = f.read()
lst = file.split('\n')[:-1]
lst2 = pas.split('\n')[:-1]
# output
out = open(pass0 + '_filteredMut','w')
whyout_output = open(whyout,'w')
dic = {}
failed = []
for i in range(len(lst)):
    info = lst[i].split('\t')
    chrom = info[0]
    pos = info[1]
    ref = lst2[i].split('\t')[3]
    alt = lst2[i].split('\t')[4]
    tRef = int(info[2])
    tAlt = int(info[3])
    nRef = int(info[4])
    nAlt = int(info[5])
    total_tumor_depth = tRef + tAlt # total tumor read depth
    vaf = float(tAlt) / total_tumor_depth # VAF for tumor
    vafbeforeout.write(str(chrom)+'\t'+str(pos)+'\t'+str(vaf)+'\t'+str(ref)+"\t"+str(alt)+'\n')
    if total_tumor_depth < 10:
        failed.append(i)
        whyout_output.write(str(chrom)+'_'+str(pos)+'\t'+ "total_tumor_depth < 10 and the death is ," + str(total_tumor_depth)+"\n")
        
    else:
        if vaf < 0.05:
            failed.append(i)
            whyout_output.write(str(chrom)+'_'+str(pos)+'\t'+"VAF < 0.05 and the VAF is ,"+ str(vaf)+ '\n')

        else:
            if tAlt <= 5 and vaf < 0.15:
                failed.append(i)
                whyout_output.write(str(chrom)+'_'+str(pos)+'\t'+"tAlt <= 5 and vaf < 0.15, and tAlt and vaf is "+ str(tAlt)+','+str(vaf)+"\n")
            else:
                if total_tumor_depth < 20 and vaf < 0.2:
                    failed.append(i)
                    whyout_output.write(str(chrom)+'_'+str(pos)+'\t'+"total_tumor_depth < 20 and vaf < 0.2 and total_tumor_depth and vaf is "+ str(total_tumor_depth)+","+str(vaf)+"\n")
                else:
                    if nAlt >= 3:
                        failed.append(i)
                        whyout_output.write(str(chrom)+'_'+str(pos)+'\t'+"nAlt >= 3 and nAlt is "+ str(nAlt)+'\n')
                    else:
                        if chrom not in dic.keys():
                            dic[chrom] = []
                        dic[chrom].append(pos)
                            
for i in range(len(lst2)):
    info = lst2[i].split('\t')
    chrom = info[0]
    pos = info[1]
    info = lst[i].split('\t')
    ref = lst2[i].split('\t')[3]
    alt = lst2[i].split('\t')[4]
    tRef = int(info[2])
    tAlt = int(info[3])
    nRef = int(info[4])
    nAlt = int(info[5])
    total_tumor_depth = tRef + tAlt # total tumor read depth
    vaf = float(tAlt) / total_tumor_depth # VAF for tumor
    if chrom in dic.keys():
        if pos in dic[chrom]:
            string = lst2[i] + '\n'
            vafafterout.write(str(chrom)+'\t'+str(pos)+'\t'+str(vaf)+'\t'+str(ref)+"\t"+str(alt)+'\n')
            out.write(string)
out.close()
vafbeforeout.close()
vafafterout.close()
whyout_output.close()
