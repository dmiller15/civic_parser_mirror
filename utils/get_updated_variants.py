import os
import sys

var_ori = {}
f = open(sys.argv[1])
header = f.readline()
for line in f:
	tmp = line.split("\t")
	var_ori[tmp[0]] = True

print header,
var_fix = {}
f = open(sys.argv[2])
header = f.readline()
for line in f:
	tmp = line.split("\t")
	if tmp[0] not in var_ori:
		print line,