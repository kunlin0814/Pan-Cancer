#!/usr/bin/python

import sys

start = sys.argv[1]
end = sys.argv[2]

string = 'qdel '

for i in range(int(start),int(end)+1):
	sample = str(i)+" "
	string += sample

print string

