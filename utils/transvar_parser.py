import os
import sys
import re

f = open(sys.argv[2], "w")
f.write("\t".join(["civic_var_id", "parse_note"]))
f.write("\n")

for line in open(sys.argv[1]):
	t = line.strip().split("\t")
	t1, t2 = t[4], t[6]
	variant["hgvs.p2g"].append(t1.split("/")[0])
	for t3 in t2.split(";"):
		t4 = t3.split("=")
		if t4[0] == "left_align_gDNA":
			variant["hgvs.p2g"].append(t4[1])
		elif t4[0] == "unalign_gDNA":
			variant["hgvs.p2g"].append(t4[1])
	f.write("\n")
