import os
import sys
import hgvs.parser

civic_var = {}
f1 = open("../data/civic_variants_extraction.tsv")
f1.readline()
for line in f1:
	tmp = line.strip().split("\t")
	civic_var[tmp[0]] = tmp[1:]

f2 = open("../data/gdc_parsed_civic_variants_all.tsv")
header = f2.readline().split("\t")
new_header = header[0:1] + ["gen_id"] + header[7:10] + header[17:20] + header[21:22]

f3 = open(sys.argv[1], "w")
f3.write("\t".join(new_header) + "\n")
hp = hgvs.parser.Parser()
for line in f2:
	tmp = line.split("\t")
	var_info = civic_var[tmp[0]]
	if tmp[17] != "" or tmp[19] != "":
		aa1_p = tmp[17]
		if tmp[17] != "":
			t = []
			for t1 in tmp[17].split(";"):
				var_p = hp.parse_hgvs_variant(t1)
				t.append(str(var_p.format(conf={"p_3_letter": False})))
			aa1_p = ";".join(t)
		f3.write("\t".join(tmp[0:1] + var_info[0:1] + tmp[7:10] + [str(aa1_p)] + tmp[18:20] + tmp[21:22]) + "\n")
f3.close()
