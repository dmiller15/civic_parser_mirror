import os
import sys
import argparse

parser = argparse.ArgumentParser(description='Generate civic variants in vcf format.')
parser.add_argument('-i', '--civic_var', type=str, required=True, help='input civic variants with cDNA info')
parser.add_argument('-t', '--transvar_var', type=str, required=True, help='transvar annotation of cDNA with vcf info')
parser.add_argument('-o', '--out_vcf', type=str, required=True, help='output variants in vcf format')
args = vars(parser.parse_args())
civic_fn = args['civic_var']
transvar_fn = args['transvar_var']
output_fn = args['out_vcf']

#cDNA_transvar_parsed_civic_variants_gseq.tsv, prediction of gDNA from cDNA by TransVar
gseq = {}
f = open(transvar_fn)
f.readline()
for line in f:
        tmp = line.strip().split("\t")
        if tmp[1] == ".":
                print("Not a proper variant:", line)
                continue
        gseq[tmp[0]] = tmp

fout = open(output_fn, "w")
fout.write("##fileformat=VCFv4.1\n")
fout.write("##fileDate=20190803\n")
fout.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

#cDNA_parsed_civic_variants.tsv, parsed cDNA from CIVIC database (which have no gDNA)
fin = open(civic_fn)
fin.readline()
for line in fin:
        cvar = line.strip().split("\t")
        if cvar[6] != "GRCh37":
                print("Not GRCh37 build or not proper contig name:", line)
                continue
        info_arr = []
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
                info_arr.append("gene_id=" + cvar[7])
        else:
                info_arr.append("gene_id=None")
        if cvar[8] != "":
                info_arr.append("entrez_id=" + cvar[8])
        else:
                info_arr.append("entrez_id=None")
        if cvar[9] != "":
                info_arr.append("entrez_name=" + cvar[9])
        else:
                info_arr.append("entrez_name=None")
        if cvar[10] != "":
                info_arr.append("civic_var_name=" + cvar[10].replace(' ', '_'))
        else:
                info_arr.append("civic_var_name=None")
        if cvar[13] != "":
                info_arr.append("hgvs.g.parsed=" + cvar[13])
        else:
                info_arr.append("hgvs.g.parsed=None")
        if cvar[15] != "":
                info_arr.append("hgvs.c.parsed=" + cvar[15])
        else:
                info_arr.append("hgvs.c.parsed=None")
        if cvar[17] != "":
                info_arr.append("hgvs.c2g.parsed=" + cvar[17])
        else:
                info_arr.append("hgvs.c2g.parsed=None")
        if cvar[19] != "":
                info_arr.append("hgvs.p.parsed=" + cvar[19])
        else:
                info_arr.append("hgvs.p.parsed=None")
        if cvar[21] != "":
                info_arr.append("vname.hgvs.p=" + cvar[21])
        else:
                info_arr.append("vname.hgvs.p=None")
        if cvar[22] != "":
                info_arr.append("vname.hgvs.p.parsed=" + cvar[22])
        else:
                info_arr.append("vname.hgvs.p.parsed=None")

        for cDNA in cvar[15].split(";"):
                if cDNA not in gseq:
                        continue
                cDNA_var = gseq[cDNA]
                info_arr.append("transvar_cDNA_transcript=" + cDNA_var[1].replace(' ', '_'))
                info_arr.append("transvar_cDNA_gene=" + cDNA_var[2].replace(' ', '_'))
                info_arr.append("transvar_cDNA_strand=" + cDNA_var[3].replace(' ', '_'))
                tmp = cDNA_var[4].split("/")
                info_arr.append("transvar_cDNA_gDNA=" + tmp[0].replace(' ', '_'))
                info_arr.append("transvar_cDNA_cDNA=" + tmp[1].replace(' ', '_'))
                info_arr.append("transvar_cDNA_protein=" + tmp[2].replace(' ', '_'))
                info_arr.append("transvar_cDNA_region=" + cDNA_var[5].replace(' ', '_'))
                info_arr.append(cDNA_var[6].strip())
                fout.write("\t".join(cDNA_var[7:9] + ["."] + cDNA_var[9:11] + [".", "."] + [";".join(info_arr)]) + "\n")

fout.close()
