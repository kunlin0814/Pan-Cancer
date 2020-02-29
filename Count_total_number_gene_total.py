#!/usr/bin/python

import sys

#file0 = sys.argv[1]

with open('/Users/kun-linho/Desktop/DNA-repair/lymphoma_mutect_DNArepair_somatic.txt','r') as f:
	file=f.read()

header = file.split('\n')[0]
lst = file.split('\n')[1:-1]

genes = header.split('\t')[1:]

def make_dic(genes):
	dic = {}
	for i in range(len(genes)):
		if genes[i] not in dic.keys():
			dic[genes[i]] = 0
	return dic

nonsyn = make_dic(genes)
totalMut = make_dic(genes)

for i in range(len(lst)):
	rec = lst[i].split('\t')[1:]
	for j in range(len(rec)):
		gene = genes[j]
		p1 = int(rec[j].split('/')[0])
		p2 = int(rec[j].split('/')[1])
		nonsyn[gene] += p1
		totalMut[gene] += p2

result = []
for i in range(len(genes)):
	gene = genes[i]
	string = str(nonsyn[gene]) + '/' + str(totalMut[gene])
	result.append(string)

string2 = 'Total\t' + '\t'.join(result)
print (string2)



