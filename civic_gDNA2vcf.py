import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Generate civic variants in vcf format.')
parser.add_argument('-c', '--contig_map', type=str, required=True, help='contig name mapping file')
parser.add_argument('-i', '--civic_var', type=str, required=True, help='input variants with gDNA info')
parser.add_argument('-t', '--transvar_var', type=str, required=True, help='transvar conversion with vcf info (use --gseq)')
parser.add_argument('-o', '--out_vcf', type=str, required=True, help='output variants in vcf format')
args = vars(parser.parse_args())
contig_map_fn = args['contig_map']
civic_fn = args['civic_var']
transvar_fn = args['transvar_var']
output_fn = args['out_vcf']

genid_chr = {}
f = open(contig_map_fn)
f.readline()
for line in f:
	tmp = line.strip().split("\t")
	genid_chr[tmp[2]] = tmp[3]
f.close()

gseq = {}
f = open(transvar_fn)
f.readline()
for line in f:
        tmp = line.strip().split("\t")
        if len(tmp) == 11:
                gseq[tmp[0]] = tmp[7:11]
        else:
                print("Not a usable prediction:", line)
f.close()

fout = open(output_fn, "w")
fout.write("##fileformat=VCFv4.1\n")
fout.write("##fileDate=20190803\n")
fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

fin = open(civic_fn)
fin.readline()
for line in fin:
	cvar = line.strip().split("\t")
	if cvar[6] != "None" and cvar[6] != "GRCh37":
		print("Not proper build GRCh37", cvar)
	else:
		info_arr = []
		info_arr.append("civic_var_id=" + cvar[0])
		info_arr.append("chr_start_stop_ref_alt=" + cvar[1])
		info_arr.append("transcript=" + cvar[2])
		info_arr.append("chr_start_stop_ref_alt_2=" + cvar[3])
		info_arr.append("transcript_2=" + cvar[4])
		info_arr.append("ensembl_version=" + cvar[5])
		info_arr.append("ref_build=" + cvar[6])
		info_arr.append("gene_id=" + cvar[7])
		info_arr.append("entrez_id=" + cvar[8])
		info_arr.append("entrez_name=" + cvar[9])
		info_arr.append("civic_var_name=" + cvar[10].replace(' ', '_'))
		info_arr.append("hgvs.g.parsed=" + cvar[13])
		info_arr.append("hgvs.c.parsed=" + cvar[15])
		info_arr.append("hgvs.c2g.parsed=" + cvar[17])
		info_arr.append("hgvs.p.parsed=" + cvar[19])
		info_arr.append("vname.hgvs.p=" + cvar[21])
		info_arr.append("vname.hgvs.p.parsed=" + cvar[22])
                
		for tmp in cvar[13].split(";"):
			try:
				tmp1 = tmp.split(":")
				gDNA = genid_chr[tmp1[0]] + ":" + tmp1[1]
				gvar = gseq[gDNA]
				fout.write("\t".join(gvar[0:2] + ["."] + gvar[2:4] + [".", "."]) + "\t" + ";".join(info_arr) + "\n")
			except:
				print(tmp, sys.exc_info()[0])
fin.close()
fout.close()
