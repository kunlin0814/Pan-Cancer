#!/usr/bin/python

import sys
import glob
import gzip

def withinRegion(pos, cds_start, cds_end):
	if (pos <= cds_end and cds_start <= pos):
		return 1
	elif (pos > cds_end):
		return 0
	elif (pos < cds_start):
		return -1


def binarySearch (arr, left, right, pos): 
	if right >= left: 
		mid = int((left+right)/2)
        	cds_start_value = int(arr[mid][0])
        	cds_end_value = int(arr[mid][1])
		regionIdx = withinRegion(pos, cds_start_value, cds_end_value)
        	if (regionIdx == 1):
        		return mid 
        	elif (regionIdx == -1):
        		return binarySearch(arr, left, mid-1, pos) 
        	elif (regionIdx == 0) : 
        		return binarySearch(arr, mid+1, right, pos) 
	else: 
        	return -2

def filter_vcf_lst(lst):
	result = []
	for i in range(len(lst)):
		if lst[i][0] != '#':
			info = lst[i].split('\t')
			status = info[6]
			if status == 'PASS':
				result.append(lst[i])
	return result

if __name__=='__main__':
	file_lst = sys.argv[1]
	CDS_file = sys.argv[2]
	file_paths = glob.glob(file_lst)
	with open(CDS_file,'r') as f:
		file = f.read()
	CDS = file.split('\n')[:-1]
	# create CDS dic
	Total_interval_dict = {}
	for i in range(len(CDS)):
		chrom = CDS[i].split(':')[0]
		#print CDS[i]
		start = int(CDS[i].split(':')[1].split('-')[0])
		end = int(CDS[i].split(':')[1].split('-')[1])
		if chrom not in Total_interval_dict.keys():
			Total_interval_dict[chrom]=[(start,end)]
		else:
			Total_interval_dict[chrom].append((start,end))
	for i in Total_interval_dict.keys():
		Total_interval_dict[i].sort()
	# read in files
	for file0 in file_paths:
		with gzip.open(file0,'rb') as f:
			file = f.read()
		lst = file.split('\n')[:-1]
		filtered_lst = filter_vcf_lst(lst)
		out = open(file0.split('.gz')[0] + '_canFam3.1.99_CDS','w')
		for i in range(len(filtered_lst)):
			rec = filtered_lst[i]
			info = rec.split('\t')
			chrom = info[0]
			pos = int(info[1])
			if chrom in Total_interval_dict.keys():
				chrom_cds_lst = Total_interval_dict[chrom]
				if (binarySearch(chrom_cds_lst,0, len(chrom_cds_lst)-1,pos)!=-2):
					string = rec + '\n'
					out.write(string)
		out.close()

