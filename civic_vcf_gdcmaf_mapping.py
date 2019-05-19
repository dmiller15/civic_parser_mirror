#!/usr/bin/env python3

import os
import sys
import argparse
import hgvs.parser

parser = argparse.ArgumentParser(description='Generate civic variants in vcf format.')
parser.add_argument('-i', '--gene_code', type=str, required=True, help='gene info from GDC-used gene model')
parser.add_argument('-g', '--gdna_var', type=str, required=True, help='input civic variants for gDNA')
parser.add_argument('-gv', '--gdna_vcf', type=str, required=True, help='input gDNA in VCF format after Liftover')
parser.add_argument('-c', '--cdna_var', type=str, required=True, help='input civic variants for cDNA')
parser.add_argument('-cv', '--cdna_vcf', type=str, required=True, help='input cDNA-predicted gDNA in VCF format after Liftover')
parser.add_argument('-p', '--prot_var', type=str, required=True, help='input civic variants for prot')
parser.add_argument('-o', '--out_dir', type=str, required=True, help='output directory for mapping/unmapping files')

args = vars(parser.parse_args())
gene_code_fn = args['gene_code']
gdna_var_fn = args['gdna_var']
gdna_vcf_fn = args['gdna_vcf']
cdna_var_fn = args['cdna_var']
cdna_vcf_fn = args['cdna_vcf']
prot_var_fn = args['prot_var']
out_dir = args['out_dir']

"""
Convert variants in VCF format to GDCMAF format.
Modify from @Kyle Hernandez's code in https://github.com/NCI-GDC/aliquot-maf-tools
"""
def vcf2maf_loc_allele(position, ref_allele, var_allele):
    ref_length, var_length = len(ref_allele), len(var_allele)

    # Remove any prefixed reference bps from all alleles, using "-" for simple indels
    while ref_allele and var_allele and ref_allele[0] == var_allele[0] and ref_allele != var_allele:
        ref_allele = ref_allele[1:] if len(ref_allele) > 1 else "-"
        var_allele = var_allele[1:] if len(var_allele) > 1 else "-"
        ref_length -= 1
        var_length -= 1
        position += 1

    # position variables
    start, stop = None, None

    # Handle SNPs, DNPs, TNPs, or anything larger
    if ref_length == var_length:
        start, stop = position, position + var_length - 1

    # Handle all indels, including those complex ones which contain substitutions
    elif ref_length != var_length:
        # Insertions
        if ref_length < var_length:
            start = position - 1 if ref_allele == "-" else position
            stop  = position if ref_allele == "-" else position + ref_length - 1
        # Deletions
        else:
            start, stop = position, position + ref_length - 1

    maf_var = {
      'start': start,
      'stop': stop,
      'ref_allele': ref_allele,
      'var_allele': var_allele,
    }
    return maf_var

def get_info(val):
    info = {}
    for v in val.split(";"):
        tmp1 = v.split("=")
        if len(tmp1) == 2:
            info[tmp1[0]] = tmp1[1]
    return info

#Get gene_name and gene_id from genecode info
genecode = {}
f = open(gene_code_fn) #gencode.gene.info.v22.tsv
f.readline()
for line in f:
    tmp = line.split("\t")
    genecode[tmp[1]] = tmp[0]

if not os.path.exists(out_dir):
	os.makedirs(out_dir)

#Mapping file for gDNA/cDNA civic info
fout = open(os.path.join(out_dir, "civic_gdcmaf_mapping_dna.tsv"), "w")
fout.write("civic_var_id\tcivic_gene_id\tsource\tchromosome\tstart_position\treference_allele\talternative_allele\n")

fin = open(gdna_vcf_fn) #gDNA_transvar_parsed_civic_variants_combined_lifted_over.vcf
for line in fin:
    if line[0] == "#":
        continue
    tmp = line.strip().split("\t")
    info = get_info(tmp[7])
    chrom, pos, ref, alt = tmp[0], int(tmp[1]), tmp[3], tmp[4]
    civic_var_id, civic_gene_id = info["civic_var_id"], info["gene_id"]
    maf_var = vcf2maf_loc_allele(pos, ref, alt)
    fout.write("\t".join([civic_var_id, civic_gene_id, "gDNA", chrom, str(maf_var["start"]), maf_var["ref_allele"], maf_var["var_allele"]]) + "\n")
fin.close()

fin = open(cdna_vcf_fn) #cDNA_transvar_parsed_civic_variants_combined_lifted_over.vcf
for line in fin:
    if line[0] == "#":
        continue
    tmp = line.strip().split("\t")
    info = get_info(tmp[7])
    chrom, pos, ref, alt = tmp[0], int(tmp[1]), tmp[3], tmp[4]
    civic_var_id, civic_gene_id = info["civic_var_id"], info["gene_id"]
    maf_var = vcf2maf_loc_allele(pos, ref, alt)
    fout.write("\t".join([civic_var_id, civic_gene_id, "cDNA", chrom, str(maf_var["start"]), maf_var["ref_allele"], maf_var["var_allele"]]) + "\n")
fin.close()

fout.close()

#Mapping file for prot civic info
fout = open(os.path.join(out_dir, "civic_gdcmaf_mapping_prot.tsv"), "w")
fout.write("civic_var_id\tcivic_gene_id\thugo_symbol\tgene\thgvs.p\tsource\n")

#To extract and process hgvs.p from all civic variants
hp = hgvs.parser.Parser()
fin = open(prot_var_fn)
fin.readline()
for line in fin:
    tmp = line.split("\t")
    civic_var_id = tmp[0]
    #if civic_var_id in gDNA_var:
    #    continue
    civic_gene_id = tmp[7]
    civic_entrez_name = tmp[9]
    aa1_p = ""
    if tmp[19] != "": #get hgvs p
        var_p = hp.parse_hgvs_variant(tmp[19])
        var_p_aa1 = str(var_p.format(conf={"p_3_letter": False}))
        aa1_p = var_p_aa1.split(":")[1]
    elif tmp[21] != "": #get vname p
        aa1_p = "p." + tmp[21]
    else:
        print("No p info:", tmp[19:24])
        continue
    gene = ""
    if civic_entrez_name in genecode:
        gene = genecode[civic_entrez_name]
    else:
        print("No geneconde mapping info:", civic_entrez_name)
    fout.write("\t".join([civic_var_id, civic_gene_id, civic_entrez_name, gene, aa1_p, "Protein"]) + "\n")

fout.close()
