import os
import sys

civic_var = {}
f1 = open(sys.argv[1])
f1.readline()
for line in f1:
	tmp = line.strip().split("\t")
	civic_var[tmp[0]] = tmp[1:]

f2 = open(sys.argv[2])
header = f2.readline().split("\t")
#new_header = header[0:7] + ["gen_id", "entrez_id"] + header[7:]
new_header = header[0:4] + ["gen_id", "entrez_id"] + header[4:]

f3 = open(sys.argv[3], "w")
f3.write("\t".join(new_header))
for line in f2:
	tmp = line.split("\t")
	#var_info = civic_var[tmp[7]]
	#f3.write("\t".join(tmp[0:7] + var_info[0:2] + tmp[7:]))
	var_info = civic_var[tmp[3]]
	f3.write("\t".join(tmp[0:4] + var_info[0:2] + tmp[4:]))
f3.close()