import os
import sys
import re
import subprocess

from civicpy import civic
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper

fix_names = {
	"Asn67fs"   : "N67fs",
	"Glu34Lys"  : "E34K",
	"EZH2 Y641F": "Y641F",
	"MLL-MLLT3" : "KMT2A-MLLT3",
	"BCR-ABL"   : "BCR-ABL1",
}

gdc_variant_types = {
	"sub"   : "point_mutation",
	"dup"   : "duplication",
	"del"   : "small_deletion",
	"ins"   : "insertion",
	"delins": "delins",
	"fs"    : "frameshift"
	#"splicing": "splicing"	
}

fusion_variant_types = {
      "FUSIONS" :"fusion_ns",
      "Fusions" :"fusion_ns",
      "fusions" : "fusion_ns",
      "FUSION"  : "fusion",
      "Fusion"  : "fusion",
      "fusion"  : "fusion"
}

other_variant_types = {
	"OVEREXPRESSION"         : "overexpression",
	"UNDEREXPRESSION"        : "underexpression",
	"DELETION"               : "large_deletion",
	"EXPRESSION"             : "expression",
	"AMPLIFICATION"          : "amplification",
	"MUTATION"               : "mutation_ns",
	"LOSS-OF-FUNCTION"       : "loss-of-function",
	"TRUNCATING MUTATION"    : "truncation",
	"FRAMESHIFT TRUNCATION"  : "frameshift_truncation",
	"PHOSPHORYLATION"        : "phosphorylation",
	"LOSS"                   : "loss",
	"METHYLATION"            : "methylation",
	"NUCLEAR EXPRESSION"     : "nuclear-expresssion",
	"CYTOPLASMIC EXPRESSION" : "cytoplasmic-expression",
	"FRAMESHIFT MUTATION"    : "frameshift_ns",
	"REARRANGEMENT"			 : "rearrangement",
	"COPY NUMBER VARIATION"  : "copy_number_change",
	"Splicing alteration"    : "splicing_ns"
}

#all_ids = range(950, 951)
all_ids = civic.get_all_variant_ids()
all_var = civic.get_variants_by_ids(all_ids)

f = open(sys.argv[1], "w")
f.write("\t".join(["civic_var_id", "chr_start_stop_ref_alt", "transcript", "chr_start_stop_ref_alt_2", "transcript_2", \
		"ensembl_version", "ref_build", "entrez_gene_name", "civic_var_name", "civic_var_types", "civic_hgvs_exp", \
		"hgvs.g.parsed", "hgvs.g.var_types", "hgvs.c.parsed", "hgvs.c.var_types", "hgvs.c2g.parsed", "hgvs.c2g.var_types", \
		"hgvs.p.parsed", "hgvs.p.var_types", "vname.hgvs.p", "vname.hgvs.p.parsed", "vname.hgvs.p.var_types", \
		"vname.hgvs.c", "vname.hgvs.c.parsed", "vname.hgvs.c.var_types", "vname.other_var_types", "parse_note"]))
f.write("\n")

hp = hgvs.parser.Parser()
for v in sorted(all_var, key=lambda x: x.id):
	print("processing variant", str(v.id))

	c = v.coordinates
	f.write("\t".join([str(v.id), str(c.chromosome) + "_" + str(c.start) + "_" + str(c.stop) + "_" + str(c.reference_bases) \
			+ "_" + str(c.variant_bases), str(c.representative_transcript), str(c.chromosome) + "_" + str(c.start2) + "_" + str(c.stop2), \
			str(c.representative_transcript2), str(c.ensembl_version), str(c.reference_build), str(v.entrez_name), str(v.name), \
			";".join(str(t.name) for t in v.variant_types), str(v.hgvs_expressions)]) + "\t")

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
				parsed_var_c2g = None
				try:
					hdp = hgvs.dataproviders.uta.connect()
					am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method='splign', replace_reference=True)
					parsed_var_c2g = am.c_to_g(parsed_var)
					variant["hgvs.c2g"].append(parsed_var_c2g)
					variant["hgvs.c2g.var_types"].append(parsed_var_c2g.posedit.edit.type)
				except:
					variant["parse_note"].append("unable to convert c to g for variant " + str(parsed_var))
				try:
					#print("using TransVar...", str(parsed_var))
				#if parsed_var_c2g == None:
					a = subprocess.run(["docker", "run", "-v", "/Users/namsyvo/.transvar.download/:/data", "-ti", \
						"zhouwanding/transvar:2.4.6", "transvar", "canno", "--ensembl", "--reference", "/data/hg19.fa", \
						"-i", str(parsed_var)], capture_output=True)
					r = str(a.stdout).split("\\n")[1].strip()
					t1 = r.split("\\t")[4]
					#print("c2g-EST, gDNA:", t1.split("/")[0])
					variant["hgvs.c2g"].append(t1.split("/")[0])
					t2 = r.split("\\t")[6]
					for t3 in t2.split(";"):
						t4 = t3.split("=")
						if t4[0] == "left_align_gDNA":
							#print("c2g-EST, left_align_gDNA:", t4[1])
							variant["hgvs.c2g"].append(t4[1])
				except:
					variant["parse_note"].append("unable to convert c to g with EST for variant " + str(parsed_var))

			elif parsed_var.type == "p":
				variant["hgvs.p"].append(parsed_var)
				variant["hgvs.p.var_types"].append(t)
		except:
			variant["parse_note"].append("unable to parse hgvs variant" + str(hgvs_var))

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
				pass
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
				pass

	if len(variant["vname.hgvs.p"]) == 0 and len(variant["vname.hgvs.c"]) == 0:
		variant["parse_note"].append("unable to parse civic variant name")

	if len(variant["vname.hgvs.p"]) == 0 and len(variant["vname.hgvs.c"]) and \
		len(variant["hgvs.g"]) == 0 and len(variant["hgvs.c2g"]) == 0 and len(variant["hgvs.p"]) == 0:
		tmp = var_name.strip()
		if tmp in other_variant_types:
			variant["vname.other_var_types"] = other_variant_types[tmp]
		else:
			variant["vname.other_var_types"] = "undetermined_var_types"

	f.write("\t".join([";".join(str(tmp) for tmp in variant["hgvs.g"]), ";".join(str(tmp) for tmp in variant["hgvs.g.var_types"]), \
			";".join(str(tmp) for tmp in variant["hgvs.c"]), ";".join(str(tmp) for tmp in variant["hgvs.c.var_types"]), \
			";".join(str(tmp) for tmp in variant["hgvs.c2g"]), ";".join(str(tmp) for tmp in variant["hgvs.c2g.var_types"]), \
			";".join(str(tmp) for tmp in variant["hgvs.p"]), ";".join(str(tmp) for tmp in variant["hgvs.p.var_types"]), \
			";".join(str(tmp) for tmp in variant["vname.hgvs.p"]), ";".join(str(tmp) for tmp in variant["vname.hgvs.p.parsed"]), \
			";".join(str(tmp) for tmp in variant["vname.hgvs.p.var_types"]), ";".join(str(tmp) for tmp in variant["vname.hgvs.c"]), \
			";".join(str(tmp) for tmp in variant["vname.hgvs.c.parsed"]), ";".join(str(tmp) for tmp in variant["vname.hgvs.p.var_types"]), \
			str(variant["vname.other_var_types"]), ";".join(str(tmp) for tmp in variant["parse_note"])]))
	f.write("\n")
