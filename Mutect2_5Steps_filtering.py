# -*- coding: utf-8 -*-
"""
Created on Mon Jun  8 20:05:36 2020

@author: abc73_000
"""

import sys
### create a Dict of Alt allele frequency from tumor samples
### Before five steps filtering


file0 = sys.argv[1]
#'G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\Mutect1_Test\\PON_DbSNP_filtered-BenGal_MuTect2_GATK4.vcf'

#'sys.argv[1]'
pass0 = sys.argv[2]
#'G:\\MAC_Research_Data\\Pan_cancer\\Pan_cancer-analysis\\Mutation_rate\\Mutect1_Test\\PON_DbSNP_filtered-BenGal_MuTect2_GATK4.vcf'
#sys.argv[2]


with open(file0,'r') as f:
	file = f.read()
        
with open(pass0,'r') as f:
	pas = f.read()
	
lst = file.split('\n')[:-1]
lst2 = pas.split('\n')[:-1]

### check normal tumor columns

def findColumn(file):
    normal_index=0
    tumor_index=0
    with open(file,'r')as f:
        for line in f:
            if (line[0]=='#' and line[1]!='#'):
                header = line.split('\t')
                for i in range(len(header)):
                    if 'normal' in header[i].lower():
                        normal_index = int(i)
                for j in range(len(header)):
                    if 'tumor' in header[j].lower():
                        tumor_index = int(j)
    return (normal_index,tumor_index)


normal_col,tumor_col = findColumn(file0)
normal_col = int(normal_col)
tumor_col = int(tumor_col)

print("normal column is ", normal_col, "Tumor column is ", tumor_col)
# output

out = open(pass0 + '_filteredMut','w')

VAF_before= sys.argv[3]
VAF_after = sys.argv[4]


VAF_out_before = open(VAF_before,'w')
VAF_out_after = open(VAF_after,'w')


dic = {}
failed = []
for i in range(len(lst)):
    if (lst[i][0] == '#'):
        pass
    else:
        if lst[i].split('\t')[6].upper() == "PASS":
            normal_info = lst[i].split('\t')[normal_col]
            tumor_info = lst[i].split('\t')[tumor_col]
            chrom = lst[i].split('\t')[0]
            pos = lst[i].split('\t')[1]
            tRef = int(tumor_info.split(':')[1].split(',')[0])
            tAlt = int(tumor_info.split(':')[1].split(',')[1])
            nRef = int(normal_info.split(':')[1].split(',')[0])
            nAlt = int(normal_info.split(':')[1].split(',')[1])
            total_tumor_depth = tRef + tAlt # total tumor read depth
            #vaf = float(tumor_info.split(":")[2])
            vaf = float(tAlt) / total_tumor_depth # VAF for tumor
            VAF_out_before.write(str(chrom)+'\t'+str(pos)+'\t'+str(vaf)+'\n')
        
            if total_tumor_depth < 10:
                failed.append(i)
            else:
                if vaf < 0.05:
                    failed.append(i)
                else:
                    if tAlt <= 5 and vaf < 0.15:
                        failed.append(i)
                    else:
                        if total_tumor_depth < 20 and vaf < 0.2:
                            failed.append(i)
                        else:
                            if nAlt >= 3:
                                failed.append(i)
                            else:
                                if chrom not in dic.keys():
                                    dic[chrom] = []
                                dic[chrom].append(pos)
                            

for i in range(len(lst2)):
    if (lst2[i][0] == '#'):
        pass
    else:
        if lst2[i].split('\t')[6].upper() == "PASS":
            info = lst2[i].split('\t')
            chrom = info[0]
            pos = info[1]
            normal_info = info[normal_col]
            tumor_info = info[tumor_col]
            tRef = int(tumor_info.split(':')[1].split(',')[0])
            tAlt = int(tumor_info.split(':')[1].split(',')[1])
            nRef = int(normal_info.split(':')[1].split(',')[0])
            nAlt = int(normal_info.split(':')[1].split(',')[1])
            total_tumor_depth = tRef + tAlt # total tumor read depth
            #vaf = float(tumor_info.split(":")[2])
            vaf = float(tAlt) / total_tumor_depth
            if chrom in dic.keys():
                if pos in dic[chrom]:
                    string = lst2[i] + '\n'
                    out.write(string)
                    VAF_out_after.write(str(chrom)+'\t'+str(pos)+'\t'+str(vaf)+'\n')
            
out.close()
VAF_out_after.close()
VAF_out_before.close()

    
    
#test = after5steps[2].split('\t')[9]
#re.search(r'(\d+:\d+,\d+:)(.:\d+:)(.)+',test).group(0)
 