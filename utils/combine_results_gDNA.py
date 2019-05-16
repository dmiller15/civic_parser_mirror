'''
Compare civic parsing (gDNA) with gdcmaf files
'''

import os
import sys

def get_info(val):
    info = {}
    for v in val.split(";"):
        tmp1 = v.split("=")
        if len(tmp1) == 2:
            info[tmp1[0]] = tmp1[1]
    return info

var_info = {}
for line in open("/data/results/cmp_civic_transvar_gdcmaf_TCGA_gDNA.tsv"):
    tmp = line.strip().split("\t")
    var_info[tmp[0]] = tmp[1]

f = open(sys.argv[3], "w")
f.write("parsed_var_hg38\tparsed_type\tvar_freq_in_gdc_maf\tcivic_var_id\thgvs.g\n")
for line in open("data/gDNA/gDNA_parsed_civic_transvar_combined_lifted_over.vcf"):
    if line[0] == "#":
        continue
    tmp = line.strip().split("\t")
    info = get_info(tmp[7])
    var = ":".join(tmp[0:2] + tmp[3:5])
    f.write("\t".join([var, "g", var_info[var], info["civic_var_id"], info["hgvs.g.parsed"]]) + "\n")

f.close()
