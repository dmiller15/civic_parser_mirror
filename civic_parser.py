#!/usr/bin/env python3

import os
import sys
import re
import json
import argparse

from civicpy import civic
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper

parser = argparse.ArgumentParser(description='Retrieve civic variant info from CIViC database')
parser.add_argument('-r', '--ref_dir', type=str, required=True, help='directory for reference files')
parser.add_argument('-o', '--out_dir', type=str, required=True, help='output directory for parsed files')

args = vars(parser.parse_args())
ref_dir = args['ref_dir']
results_dir = args['out_dir']

if not os.path.exists(ref_dir):
        print("Please indicate the reference directory")
        sys.exit(1)
if not os.path.exists(results_dir):
	os.makedirs(results_dir)

with open(os.path.join(ref_dir, "fix_names.json")) as f:
	fix_names = json.load(f)
with open(os.path.join(ref_dir, "gdc_variant_types.json")) as f:
	gdc_variant_types = json.load(f)
with open(os.path.join(ref_dir, "fusion_variant_types.json")) as f:
	fusion_variant_types = json.load(f)
with open(os.path.join(ref_dir, "other_variant_types.json")) as f:
	other_variant_types = json.load(f)

genid_chr = {}
f = open(os.path.join(ref_dir, "geneid_conversion.tsv"))
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	genid_chr[tmp[2]] = tmp[3]

all_ids = civic.get_all_variant_ids()
all_var = civic.get_variants_by_ids(all_ids)

f_a = open(os.path.join(results_dir, "all_parsed_civic_variants.tsv"), "w")
f_g = open(os.path.join(results_dir, "gDNA_parsed_civic_variants.tsv"), "w")
f_c = open(os.path.join(results_dir, "cDNA_parsed_civic_variants.tsv"), "w")
f_p = open(os.path.join(results_dir, "prot_parsed_civic_variants.tsv"), "w")
f_o = open(os.path.join(results_dir, "gdc_unsupported_civic_variants.tsv"), "w")
f_u = open(os.path.join(results_dir, "unparsed_civic_variants.tsv"), "w")
f_ct = open(os.path.join(results_dir, "cDNA_parsed_civic_variants_clist.tsv"), "w")
f_gt = open(os.path.join(results_dir, "gDNA_parsed_civic_variants_glist.tsv"), "w")

header = "\t".join(["civic_var_id", "chr_start_stop_ref_alt", "transcript", "chr_start_stop_ref_alt_2", "transcript_2", \
		"ensembl_version", "ref_build", "gene_id", "entrez_id", "entrez_name", "civic_var_name", "civic_var_types", "civic_hgvs_exp", \
		"hgvs.g.parsed", "hgvs.g.var_types", "hgvs.c.parsed", "hgvs.c.var_types", "hgvs.c2g.parsed", "hgvs.c2g.var_types", \
		"hgvs.p.parsed", "hgvs.p.var_types", "vname.hgvs.p", "vname.hgvs.p.parsed", "vname.hgvs.p.var_types", \
		"vname.hgvs.c", "vname.hgvs.c.parsed", "vname.hgvs.c.var_types", "vname.other_var_types", "parse_note"])
f_a.write(header + "\n")
f_g.write(header + "\n")
f_c.write(header + "\n")
f_p.write(header + "\n")
f_o.write(header + "\n")
f_u.write(header + "\n")

hp = hgvs.parser.Parser()
for v in sorted(all_var, key=lambda x: x.id):
	print("processing variant", str(v.id))

	variant = {}

	variant["hgvs.g"] = []
	variant["hgvs.c"] = []
	variant["hgvs.c2g"] = []
	variant["hgvs.p"] = []

	variant["hgvs.g.var_types"] = []
	variant["hgvs.c.var_types"] = []
	variant["hgvs.c2g.var_types"] = []
	variant["hgvs.p.var_types"] = []

	variant["vname.hgvs.p"] = []
	variant["vname.hgvs.p.parsed"] = []
	variant["vname.hgvs.p.var_types"] = []

	variant["vname.hgvs.c"] = []
	variant["vname.hgvs.c.parsed"] = []
	variant["vname.hgvs.c.var_types"] = []

	variant["vname.other_var_types"] = ""
	variant["parse_note"] = []

	for hgvs_var in v.hgvs_expressions:
		try:
			parsed_var = hp.parse_hgvs_variant(hgvs_var)
			t = parsed_var.posedit.edit.type
			if parsed_var.type == "g":
				variant["hgvs.g"].append(parsed_var)
				variant["hgvs.g.var_types"].append(t)
			elif parsed_var.type == "c":
				variant["hgvs.c"].append(parsed_var)
				variant["hgvs.c.var_types"].append(t)
				'''
				parsed_var_c2g = None
				try:
					hdp = hgvs.dataproviders.uta.connect()
					am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method='splign', replace_reference=True)
					parsed_var_c2g = am.c_to_g(parsed_var)
					variant["hgvs.c2g"].append(parsed_var_c2g)
					variant["hgvs.c2g.var_types"].append(parsed_var_c2g.posedit.edit.type)
				except:
					variant["parse_note"].append("unable to convert c to g for variant " + str(parsed_var))
					print("c2g conversion error", sys.exc_info()[0])
				'''
			elif parsed_var.type == "p":
				variant["hgvs.p"].append(parsed_var)
				variant["hgvs.p.var_types"].append(t)
			else:
				print("Unprocessed parsed_var.type: ", parsed_var.type)
		except:
			variant["parse_note"].append("unable to parse hgvs variant" + str(hgvs_var))
			print("hgvs parsing error", sys.exc_info()[0])

	var_name = str(v.name)
	vname_arr = re.split('\s|\(|\)|\+', var_name.strip())
	for vname in vname_arr:
			parsed_vname = re.sub("^p\\.", "", vname) # remove p. from begining because CIVIC is not consistent on this
			parsed_vname = re.sub("^([ACDEFGHIKLMNPQRSTVWY\\*]\\d+)X$", "\\1*", parsed_vname) # change terminal X into *
			parsed_vname = parsed_vname.replace("Ter", "*") # change Ter into *
			parsed_vname = parsed_vname.replace("FS", "fs") # change FS into fs
			parsed_vname = parsed_vname.replace("DUP", "dup") # change DUP into dup
			parsed_vname = parsed_vname.replace("DEL", "del") # change DEL into del
			parsed_vname = parsed_vname.replace("INS", "ins") # change INS into ins
			parsed_vname = parsed_vname.replace("DELins", "delins") # change DELins into delins
			
			if parsed_vname in fix_names:
				parsed_vname = fix_names[parsed_vname]
				variant["parse_note"].append("fix civic var name for vname.hgvs.p")
			if str(v.entrez_name) == "JAK2" and parsed_vname == "F694L":
				parsed_vname = "F694L"
				variant["parse_note"].append("fix civic var name for vname.hgvs.p")
			try:
				pv = hp.parse_p_posedit(parsed_vname)
				variant["vname.hgvs.p"].append(str(parsed_vname))
				variant["vname.hgvs.p.parsed"].append(str(pv))
				t = pv.edit.type
				if t in gdc_variant_types:
					variant["vname.hgvs.p.var_types"].append(gdc_variant_types[t])
				else:
					variant["vname.hgvs.p.var_types"].append(t)
			except:
				print("parse_p_posedit error", parsed_vname, sys.exc_info()[0])
			try:
				cv = hp.parse_c_posedit(parsed_vname)
				variant["vname.hgvs.c"].append(str(parsed_vname))
				variant["vname.hgvs.c.parsed"].append(str(cv))
				t = cv.edit.type
				if t in gdc_variant_types:
					variant["vname.hgvs.c.var_types"].append(gdc_variant_types[t])
				else:
					variant["vname.hgvs.c.var_types"].append(t)
			except:
				print("parse_c_posedit error", parsed_vname, sys.exc_info()[0])

	if len(variant["vname.hgvs.p"]) == 0 and len(variant["vname.hgvs.c"]) == 0:
		variant["parse_note"].append("unable to parse civic variant name")

	if len(variant["vname.hgvs.p"]) == 0 and len(variant["vname.hgvs.c"]) and \
		len(variant["hgvs.g"]) == 0 and len(variant["hgvs.c2g"]) == 0 and len(variant["hgvs.p"]) == 0:
		tmp = var_name.strip()
		if tmp in other_variant_types:
			variant["vname.other_var_types"] = other_variant_types[tmp]
		else:
			variant["vname.other_var_types"] = "undetermined_var_types"

	c = v.coordinates
	datum1 = "\t".join([str(v.id), str(c.chromosome) + "_" + str(c.start) + "_" + str(c.stop) + "_" + str(c.reference_bases) \
			+ "_" + str(c.variant_bases), str(c.representative_transcript), str(c.chromosome) + "_" + str(c.start2) + "_" + str(c.stop2), \
			str(c.representative_transcript2), str(c.ensembl_version), str(c.reference_build), str(v.gene_id), str(v.entrez_id), \
			str(v.entrez_name), str(v.name), ";".join(str(t.name) for t in v.variant_types), str(v.hgvs_expressions)]) + "\t"

	datum2 = "\t".join([";".join(str(tmp) for tmp in variant["hgvs.g"]), ";".join(str(tmp) for tmp in variant["hgvs.g.var_types"]), \
			";".join(str(tmp) for tmp in variant["hgvs.c"]), ";".join(str(tmp) for tmp in variant["hgvs.c.var_types"]), \
			";".join(str(tmp) for tmp in variant["hgvs.c2g"]), ";".join(str(tmp) for tmp in variant["hgvs.c2g.var_types"]), \
			";".join(str(tmp) for tmp in variant["hgvs.p"]), ";".join(str(tmp) for tmp in variant["hgvs.p.var_types"]), \
			";".join(str(tmp) for tmp in variant["vname.hgvs.p"]), ";".join(str(tmp) for tmp in variant["vname.hgvs.p.parsed"]), \
			";".join(str(tmp) for tmp in variant["vname.hgvs.p.var_types"]), ";".join(str(tmp) for tmp in variant["vname.hgvs.c"]), \
			";".join(str(tmp) for tmp in variant["vname.hgvs.c.parsed"]), ";".join(str(tmp) for tmp in variant["vname.hgvs.p.var_types"]), \
			str(variant["vname.other_var_types"]), ";".join(str(tmp) for tmp in variant["parse_note"])])
	f_a.write(datum1 + datum2 + "\n")
	vtypes = [str(t.name) for t in v.variant_types]
	if len(variant["hgvs.g"]) > 0:
		if "transcript_fusion" in vtypes:
			f_o.write(datum1 + datum2 + "\n")
		else:
			f_g.write(datum1 + datum2 + "\n")
			for g in variant["hgvs.g"]:
				tmp = str(g).split(":")
				if tmp[0] in genid_chr:
					g_chr = genid_chr[tmp[0]] + ":" + tmp[1]
					f_gt.write(g_chr + "\n")
	elif len(variant["hgvs.c"]) > 0:
		if "transcript_fusion" in vtypes:
			f_o.write(datum1 + datum2 + "\n")
		else:
			f_c.write(datum1 + datum2 + "\n")
			for c in variant["hgvs.c"]:
				f_ct.write(str(c) + "\n")
	elif len(variant["hgvs.p"]) > 0 or len(variant["vname.hgvs.p"]) > 0:
		if "transcript_fusion" in vtypes:
			f_o.write(datum1 + datum2 + "\n")
		else:
			f_p.write(datum1 + datum2 + "\n")
	else:
		f_u.write(datum1 + datum2 + "\n")

f_a.close()
f_g.close()
f_c.close()
f_p.close()
f_o.close()
f_u.close()
f_ct.close()
f_gt.close()
