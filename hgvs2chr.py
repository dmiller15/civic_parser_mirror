import os
import sys

genid_chr = {}
f = open(sys.argv[1])
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	genid_chr[tmp[2]] = tmp[3]

f = open(sys.argv[2])
f.readline()
c = 0
for line in f:
	tmp = line.strip().split("\t")
	if tmp[11] != "":
		c += 1
		tmp2 = tmp[11].split(":")
		try:
			print(genid_chr[tmp2[0]] + ":" + tmp2[1])
		except:
			print("error", sys.exc_info()[0])
