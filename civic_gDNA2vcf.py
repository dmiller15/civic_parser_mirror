import os
import sys

genid_chr = {}
f = open("data/geneid_conversion.tsv")
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	genid_chr[tmp[2]] = tmp[3]

civic_var = {}
f = open("data/gDNA/gDNA_parsed_civic_variants.tsv")
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	tmp2 = tmp[11].split(":")
	try:
		gDNA = genid_chr[tmp2[0]] + ":" + tmp2[1]
		civic_var[gDNA] = tmp
	except:
		print tmp2, sys.exc_info()[0]

fout = open(sys.argv[1], "w")
fout.write("##fileformat=VCFv4.1\n")
fout.write("##fileDate=20190803\n")
fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

f = open("data/gDNA/gDNA_gseq_parsed_civic_variants.tsv")
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	try:
		info_arr = []
		cvar = civic_var[tmp[0]]
		info_arr.append("civic_var_id=" + cvar[0])
		info_arr.append("chr_start_stop_ref_alt=" + cvar[1])
		info_arr.append("transcript=" + cvar[2])
		info_arr.append("chr_start_stop_ref_alt_2=" + cvar[3])
		info_arr.append("transcript_2=" + cvar[4])
		info_arr.append("ensembl_version=" + cvar[5])
		info_arr.append("ref_build=" + cvar[6])
		info_arr.append("entrez_gene_name=" + cvar[7])
		info_arr.append("civic_var_name=" + cvar[8].replace(' ', '_'))
		info_arr.append("hgvs.g.parsed=" + cvar[11])
		info_arr.append("hgvs.c.parsed=" + cvar[13])
		info_arr.append("hgvs.c2g.parsed=" + cvar[15])
		info_arr.append("hgvs.p.parsed=" + cvar[17])
		info_arr.append("vname.hgvs.p=" + cvar[19])
		info_arr.append("vname.hgvs.p.parsed=" + cvar[20])

		info_arr.append(tmp[6])
		info_arr.append("region=" + tmp[5])

		fout.write("\t".join(tmp[7:9] + ["."] + tmp[9:11] + [".", "."]) + "\t" + ";".join(info_arr) + "\n")
	except:
		print tmp[0], sys.exc_info()[0]

fout.close()
