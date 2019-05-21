import os
import sys

civic_var = {}
f1 = open(sys.argv[1]) #data/results/civic_variants_extraction.tsv
f1.readline()
for line in f1:
	tmp = line.strip().split("\t")
	civic_var[tmp[0]] = tmp[1:]

f2 = open(sys.argv[2]) #data/results2/mapping/civic_gdcmaf_mapping_dna.tsv
header = f2.readline()

f3 = open(sys.argv[3], "w") #data/results2/mapping/civic_gdcmaf_mapping_dna_updated.tsv
f3.write(header)
for line in f2:
	tmp = line.split("\t")
	var_info = civic_var[tmp[0]]
	f3.write("\t".join(tmp[0:1] + var_info[0:1] + tmp[2:]))
f3.close()
