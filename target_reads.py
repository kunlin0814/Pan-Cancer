#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 29 13:38:09 2020

@author: kun-linho
"""
from collections import Counter

# for each chromosome
# min = 23840680 
# max = 123780325
#import sys
## use 50 million bp as a cut off


def binarySearch (arr, left, right, reads_position): 
  
    # Check base case 
    if right >= left: 
  
        mid = int(left + (right - left)/2)
  
        # If element is present at the middle itself 
        if arr[mid] == reads_position: 
            return mid 
          
        # If element is smaller than mid, then it can only 
        # be present in left subarray 
        elif arr[mid] > reads_position: 
            return binarySearch(arr, left, mid-1, reads_position) 
  
        # Else the element can only be present in right subarray 
        else: 
            return binarySearch(arr, mid+1, right, reads_position) 
  
    else: 
        # Element is not present in the array 
        return -1

read_length = 101
small_interval_dict ={}
large_interval_dict ={}
interval_size ={}
with open('/Volumes/Research_Data/Pan_cancer/Mapping_source/uniq_Canis_familiaris.CanFam3.1.81.gtf-chr1-38X-CDS-forDepthOfCoverage.interval_list', 'r') as f:
    file = f.read()


CDS = file.split('\n')[:-1]
    
for i in range(len(CDS)):
    chrom = CDS[i].split(':')[0]
    start = int(CDS[i].split(':')[1].split('-')[0])
    end = int(CDS[i].split(':')[1].split('-')[1])
    if (start >= 50000000 and end >= 50000000):
        if chrom not in large_interval_dict.keys():
            large_interval_dict[chrom]=[(start, end)]
        else:
            large_interval_dict[chrom].append((start, end))
    else: 
         if chrom not in small_interval_dict.keys():
            small_interval_dict[chrom]=[(start, end)]
         else:
            small_interval_dict[chrom].append((start, end))

for i in small_interval_dict.keys():
    small_interval_dict[i].sort()
    large_interval_dict[i].sort()


        
#for i in small_interval_dict.keys():
#    chrom_end = max(small_interval_dict[i])[1]
#    chrom_start = min(small_interval_dict[i])[0]
#    distance = chrom_end-chrom_start
#    interval_size[i]= distance

#sam_file=sys.argv[1]
#output_name=sys.argv[2]
## criteria is Reads start site <=Exon end site and exon start site <= reads end site
    
unique = 0 #3
duplicate = 0 #3
Onemapped = 0 #5,9
incorrect = 0 #1
unmapped = 0 #13
#with open('C:/Users/abc73_000/Desktop/test.sam','r') as f1:
#    file = f.read()
#    sam_file = file.split('\n')[40:-1]
    
#reads_name = sam_file[0].split('\t')[0]    
#reads_chr= sam_file[0].split('\t')[2]
#reads_position = int(sam_file[0].split('\t')[3])
## criteria is Reads start site <=Exon end site and exon start site <= reads end site
transcript_list =[]    
total = 0 
pass_line =0
with open('C:/Users/abc73_000/Desktop/test.sam','r') as f1:   
    for line in f1:
        file_lst = line.split('\t')
        if '@' in file_lst[0]:
            pass
        else :
            total += 1
            reads_name = file_lst[0]
            reads_chr = file_lst[2]
            reads_position = int(file_lst[3])
            status = int(file_lst[1])%16
            if status == 3:
                for ele in file_lst:
                    if 'XT:' in ele:
                        status2 = ele.split(':')[2]
                        if status2 == 'U' or status2 == 'M':
                            unique += 1
                            if reads_position < 50000000:
                                exome_loc = small_interval_dict[reads_chr]
                                
                                for i in exome_loc :
                                    if (reads_position<= i[1] and i[0]<=(reads_position+read_length)):
                                        transcript_list.append(reads_name)
                                    
                            else:
                                 exome_loc = large_interval_dict[reads_chr]
                                 
                                 for i in exome_loc :
                                    if (reads_position<= i[1] and i[0]<=(reads_position+read_length)):
                                        transcript_list.append(reads_name)

                        
counter = Counter(transcript_list)
dupes = [key for (key, value) in counter.items() if value > 1 and key]

# total = unique + duplicate + Onemapped + incorrect + unmapped
pairs = total / 2
print ('total reads pair are ' + str(pairs) )
print ('the duplicate transcript is '+ str(dupes))
print ('uniq coordinate correctly mapped reads is ' + str(unique))
       

    
