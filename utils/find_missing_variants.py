import os
import sys

missing_var = {}
f = open(sys.argv[1])
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	missing_var[tmp[7]] = tmp[6]

f = open(sys.argv[2])
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	if tmp[0] in missing_var:
		print line.strip()

'''
gseq_var = {}
f = open(sys.argv[1])
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	tmp1 = tmp[0].split(":")[1]
	gseq_var[tmp1] = True

f = open(sys.argv[2])
print f.readline()
for line in f:
	tmp = line.strip().split("\t")
	tmp1 = tmp[11].split(":")[1]
	if tmp1 not in gseq_var:
		print line
'''