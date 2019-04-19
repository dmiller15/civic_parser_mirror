'''
Compare civic parsing (prot) with gdcmaf files
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
for line in open("/data/results/cmp_civic_transvar_gdcmaf_TCGA_prot.tsv"):
    tmp = line.strip().split("\t")
    var_info[tmp[0]] = tmp[1]

f = open(sys.argv[3], "w")
f.write("parsed_var_hg38\tparsed_type\tvar_freq_in_gdc_maf\tcivic_var_id\ttransvar_hgvs.g_hg19\tcivic_transcript\ttransvar_transcript\tcandidate_codons\taliases\tdb_sources\tisoform_match\n")
for line in open("data/prot/prot_parsed_civic_transvar_combined_lifted_over.vcf"):
    if line[0] == "#":
        continue
    tmp = line.strip().split("\t")
    info = get_info(tmp[7])
    var = ":".join(tmp[0:2] + tmp[3:5])
    tmp2 = info["transvar_prot_transcript"].split("_")
    transvar_prot_transcript = "_".join(tmp2[:-2])
    candidate_codons = info.get("candidate_codons", "None")
    aliases = info.get("aliases", "None")
    if info["transcript"] == "None" or info["transcript"] == ".":
        f.write("\t".join([var, "p", var_info[var], info["civic_var_id"], info["transvar_prot_gDNA"], info["transcript"], transvar_prot_transcript, candidate_codons, aliases, info["source"][:-1], "No_civic_info"]) + "\n")
    else:
        if transvar_prot_transcript in info["transcript"]:
            f.write("\t".join([var, "p", var_info[var], info["civic_var_id"], info["transvar_prot_gDNA"], info["transcript"], transvar_prot_transcript, candidate_codons, aliases, info["source"][:-1], "Match"]) + "\n")
        elif "ENST" not in transvar_prot_transcript:
            f.write("\t".join([var, "p", var_info[var], info["civic_var_id"], info["transvar_prot_gDNA"], info["transcript"], transvar_prot_transcript, candidate_codons, aliases, info["source"][:-1], "No_civic_info"]) + "\n")
        else:
            f.write("\t".join([var, "p", var_info[var], info["civic_var_id"], info["transvar_prot_gDNA"], info["transcript"], transvar_prot_transcript, candidate_codons, aliases, info["source"][:-1], "Not_match"]) + "\n")

f.close()
