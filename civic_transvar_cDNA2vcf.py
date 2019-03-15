import os
import sys

#cDNA_parsed_civic_variants.tsv
#parsed cDNA from CIVIC database (which have no gDNA)
civic_var = {}
f = open(sys.argv[1])
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	for cDNA_var in tmp[13].split(";"):
		civic_var[cDNA_var] = tmp

#cDNA_transvar_parsed_civic_variants.tsv
#prediction of gDNA from cDNA by TransVar
transvar_cDNA_var = {}
f = open(sys.argv[2])
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	if tmp[4] not in transvar_cDNA_var:
		transvar_cDNA_var[tmp[4]] = []
	transvar_cDNA_var[tmp[4]].append(tmp)

fout = open(sys.argv[4], "w")
fout.write("##fileformat=VCFv4.1\n")
fout.write("##fileDate=20190803\n")
fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

#cDNA_transvar_parsed_civic_gseq.tsv
#conversion of HGVS to VCF by TransVar
f = open(sys.argv[3])
f.readline()
for line in f:
	tmp = line.strip().split("\t")

	for cDNA_var in transvar_cDNA_var[tmp[0]]:
		info_arr = []
		info_arr.append("transvar_cDNA_transcript=" + cDNA_var[1].replace(' ', '_'))
		info_arr.append("transvar_cDNA_gene=" + cDNA_var[2].replace(' ', '_'))
		info_arr.append("transvar_cDNA_strand=" + cDNA_var[3].replace(' ', '_'))
		info_arr.append("transvar_cDNA_gDNA=" + cDNA_var[4].replace(' ', '_'))
		info_arr.append("transvar_cDNA_cDNA=" + cDNA_var[5].replace(' ', '_'))
		info_arr.append("transvar_cDNA_protein=" + cDNA_var[6].replace(' ', '_'))
		info_arr.append("transvar_cDNA_region=" + cDNA_var[7].replace(' ', '_'))
		info_arr.append(cDNA_var[8].strip())

		cvar = civic_var[cDNA_var[0]]
		info_arr.append("civic_var_id=" + cvar[0])
		info_arr.append("chr_start_stop_ref_alt=" + cvar[1])
		info_arr.append("transcript=" + cvar[2])
		info_arr.append("chr_start_stop_ref_alt_2=" + cvar[3])
		info_arr.append("transcript_2=" + cvar[4])

		if cvar[5] != "":
			info_arr.append("ensembl_version=" + cvar[5])
		else:
			info_arr.append("ensembl_version=None")

		if cvar[6] != "":
			info_arr.append("ref_build=" + cvar[6])
		else:
			info_arr.append("ref_build=None")

		if cvar[7] != "":
			info_arr.append("entrez_gene_name=" + cvar[7])
		else:
			info_arr.append("entrez_gene_name=None")

		if cvar[8] != "":
			info_arr.append("civic_var_name=" + cvar[8].replace(' ', '_'))
		else:
			info_arr.append("civic_var_name=None")

		if cvar[11] != "":
			info_arr.append("hgvs.g.parsed=" + cvar[11])
		else:
			info_arr.append("hgvs.g.parsed=None")

		if cvar[13] != "":
			info_arr.append("hgvs.c.parsed=" + cvar[13])
		else:
			info_arr.append("hgvs.c.parsed=None")

		if cvar[15] != "":
			info_arr.append("hgvs.c2g.parsed=" + cvar[15])
		else:
			info_arr.append("hgvs.c2g.parsed=None")

		if cvar[17] != "":
			info_arr.append("hgvs.p.parsed=" + cvar[17])
		else:
			info_arr.append("hgvs.p.parsed=None")

		if cvar[19] != "":
			info_arr.append("vname.hgvs.p=" + cvar[19])
		else:
			info_arr.append("vname.hgvs.p=None")

		if cvar[20] != "":
			info_arr.append("vname.hgvs.p.parsed=" + cvar[20])
		else:
			info_arr.append("vname.hgvs.p.parsed=None")

		info_arr.append(tmp[6])
		info_arr.append("region=" + tmp[5].strip())

		fout.write("\t".join(tmp[7:9] + ["."] + tmp[9:11] + [".", "."] + [";".join(info_arr)]) + "\n")

fout.close()
