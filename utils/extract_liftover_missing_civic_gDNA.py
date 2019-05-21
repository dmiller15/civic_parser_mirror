import os
import sys

civic_var = {}
#f = open(sys.argv[1])
f = open("data/results/cmp_civic_transvar_gdcmaf_TCGA_gDNA_full_info_update.tsv")
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	civic_var[tmp[3]] = True

#fout = open(sys.argv[3], "w")
fout = open("data/results/gDNA_parsed_civic_variants_not_liftover.tsv", "w")

#f = open(sys.argv[2])
f = open("data/gDNA/gDNA_parsed_civic_variants.tsv")
fout.write(f.readline())
for line in f:
	tmp = line.strip().split("\t")
	if tmp[0] not in civic_var:
		fout.write(line)
fout.close()
