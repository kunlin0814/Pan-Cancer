#!/usr/bin/python

import sys
import re

file=sys.argv[1]
name_pair=sys.argv[2]

with open(file,'r') as f:	
	file_info=f.read()

with open(name_pair,'r') as f:
	name_pair=f.read()

output=open(file+'_WithGeneName','w')

signature=open(file+'_Signature','w')

file_lst=file_info.split('\n')[:-1]

pair_lst=name_pair.split('\n')[:-1]

for i in range(len(file_lst)):
	info=file_lst[i].split('\t',3)
	if len(info)>=4:
		if len(info[2].split(':'))>=4:
			ensembl_id=info[2].split(':')[0]
			if len(info[2].split(':')[3].split('.'))==2:
				aa_subs=info[2].split(':')[3].split('.')[1]
				p = re.compile(r'\d+')
				sig='-'.join(p.split(aa_subs))
				signature.write(sig)
				signature.write('\n')
				for j in range(len(pair_lst)):
					ens=pair_lst[j].split('\t')[0]
					gene=pair_lst[j].split('\t')[1]
					if ensembl_id==ens:
						lst=[]
						lst.append(info[0])
						lst.append(info[1])
						lst.append(gene)
						lst.append(info[2])
						lst.append(info[3])
						output.write('\t'.join(lst))
						output.write('\n')



output.close()
signature.close()



