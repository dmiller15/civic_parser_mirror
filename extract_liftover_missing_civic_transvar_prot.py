import os
import sys

liftover_var = {}
#f = open(sys.argv[1])
f = open("data/results/cmp_civic_transvar_gdcmaf_TCGA_prot_full_info_update.tsv")
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	liftover_var[tmp[3]] = True

civic_transvar_var = {}
f = open("data/prot/prot_transvar_parsed_civic_variants.tsv")
tv_header = f.readline().strip().split("\t")[1:]
for line in f:
	tmp = line.strip().split("\t")
	if len(tmp) < 5 or "chr" not in tmp[4]:
		continue
	if tmp[0] not in civic_transvar_var:
		civic_transvar_var[tmp[0]] = []
	civic_transvar_var[tmp[0]].append(tmp[1:])

item_num = max([len(v) for _, v in civic_transvar_var.iteritems()])

#fout = open(sys.argv[3], "w")
fout = open("data/results/prot_parsed_civic_variants_liftover_missing.tsv", "w")

#f = open(sys.argv[2])
f = open("data/prot/prot_parsed_civic_variants.tsv")
header = f.readline().strip()

fout.write(header + "\thgvs.g_convertable")
for i in range(item_num):
	fout.write("\t" + "\t".join(["transvar_info_" + str(i + 1) + ":" + "|".join(tv_header)]))
fout.write("\n")

for line in f:
	tmp = line.strip().split("\t")
	if tmp[0] not in liftover_var:
		prot = tmp[7] + ":p." + tmp[19]
		if prot in civic_transvar_var:
			fout.write("\t".join(tmp) + "\tTrue")
			for v in civic_transvar_var[prot]:
				fout.write("\t" + "|".join(v))
			fout.write("\n")
		else:
			fout.write(line.strip() + "\tFalse\n")

fout.close()

