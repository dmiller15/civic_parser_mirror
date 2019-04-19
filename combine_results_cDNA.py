'''
Compare civic parsing (cDNA) with gdcmaf files
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
for line in open("/data/results/cmp_civic_transvar_gdcmaf_TCGA_cDNA.tsv"):
    tmp = line.strip().split("\t")
    var_info[tmp[0]] = tmp[1]

f = open(sys.argv[3], "w")
f.write("parsed_var_hg38\tparsed_type\tvar_freq_in_gdc_maf\tcivic_var_id\ttransvar_hgvs.g_hg19\tcivic_transcript\ttransvar_transcript\taliases\tdb_sources\tisoform_match\n")
for line in open("data/cDNA/cDNA_parsed_civic_transvar_combined_lifted_over.vcf"):
    if line[0] == "#":
        continue
    tmp = line.strip().split("\t")
    if len(tmp[3]) > 50:
        continue
    info = get_info(tmp[7])
    var = ":".join(tmp[0:2] + tmp[3:5])
    tmp2 = info["transvar_cDNA_transcript"].split("_")
    transvar_cDNA_transcript = "_".join(tmp2[:-2])
    aliases = info.get("aliases", "None")
    if info["transcript"] == "None" or info["transcript"] == ".":
        f.write("\t".join([var, "c", var_info[var], info["civic_var_id"], info["transvar_cDNA_gDNA"], info["transcript"], transvar_cDNA_transcript, aliases, info["source"], "No_civic_info"]) + "\n")
    else:
        if transvar_cDNA_transcript in info["transcript"]:
            f.write("\t".join([var, "c", var_info[var], info["civic_var_id"], info["transvar_cDNA_gDNA"], info["transcript"], transvar_cDNA_transcript, aliases, info["source"], "Match"]) + "\n")
        elif "ENST" not in transvar_cDNA_transcript or ("ENST" in info["transcript"] and len(info["transcript"]) < 15):
            f.write("\t".join([var, "c", var_info[var], info["civic_var_id"], info["transvar_cDNA_gDNA"], info["transcript"], transvar_cDNA_transcript, aliases, info["source"], "No_civic_info"]) + "\n")
        else:
            f.write("\t".join([var, "c", var_info[var], info["civic_var_id"], info["transvar_cDNA_gDNA"], info["transcript"], transvar_cDNA_transcript, aliases, info["source"], "Not_match"]) + "\n")

f.close()
