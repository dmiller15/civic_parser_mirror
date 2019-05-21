#!/usr/bin/env python3

import os
import sys
import argparse
import hgvs.parser

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

parser = argparse.ArgumentParser(description='Generate civic variant info to be used with GDC mutation indexing')
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

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

#Get gene_name and gene_id from genecode info
genecode = {}
f = open(gene_code_fn) #gencode.gene.info.v22.tsv
f.readline()
for line in f:
    tmp = line.split("\t")
    genecode[tmp[1]] = tmp[0]

#Extract gDNA from liftovered variants
long_var_id = []
gdna_vcf = {}
fin = open(gdna_vcf_fn) #gDNA_transvar_parsed_civic_variants_combined_lifted_over.vcf
for line in fin:
    if line[0] == "#":
        continue
    tmp = line.strip().split("\t")
    info = get_info(tmp[7])
    chrom, pos, ref, alt = tmp[0], int(tmp[1]), tmp[3], tmp[4]
    civic_var_id, civic_gene_id = info["civic_var_id"], info["gene_id"]
    maf_var = vcf2maf_loc_allele(pos, ref, alt)
    if len(maf_var["ref_allele"]) > 50 or len(maf_var["var_allele"]) > 50:
        long_var_id.append(civic_var_id)
    else:
        gdna_vcf[civic_var_id] = [civic_var_id, civic_gene_id, "gDNA", chrom, \
                str(maf_var["start"]), maf_var["ref_allele"], maf_var["var_allele"]]
fin.close()

#Extract cDNA from liftovered variants
cdna_vcf = {}
fin = open(cdna_vcf_fn) #cDNA_transvar_parsed_civic_variants_combined_lifted_over.vcf
for line in fin:
    if line[0] == "#":
        continue
    tmp = line.strip().split("\t")
    info = get_info(tmp[7])
    chrom, pos, ref, alt = tmp[0], int(tmp[1]), tmp[3], tmp[4]
    civic_var_id, civic_gene_id = info["civic_var_id"], info["gene_id"]
    maf_var = vcf2maf_loc_allele(pos, ref, alt)
    if len(maf_var["ref_allele"]) > 50 or len(maf_var["var_allele"]) > 50:
        long_var_id.append(civic_var_id)
    else:
        cdna_vcf[civic_var_id] = [civic_var_id, civic_gene_id, "cDNA", chrom, str(maf_var["start"]), \
                                  maf_var["ref_allele"], maf_var["var_allele"]]
fin.close()

#Extract hgvs.p from parsed civic variants
hgvsp = {}

f_umap_long_var = open(os.path.join(out_dir, "parsed_civic_variants.unmapped_long_var.tsv"), "w")

gdna_var = {}
fin = open(gdna_var_fn)
header = fin.readline()
fout = open(gdna_var_fn + ".unliftovered", "w")
fout.write(header)
f_umap_long_var.write(header)
for line in fin:
    cvar = line.split("\t")
    cvar_id = cvar[0]
    gdna_var[cvar_id] = cvar
    if cvar_id not in gdna_vcf:
        hgvsp[cvar_id] = cvar
        fout.write(line)
    if cvar_id in long_var_id:
        f_umap_long_var.write(line)
fin.close()
fout.close()

cdna_var = {}
fin = open(cdna_var_fn)
fout = open(cdna_var_fn + ".unliftovered", "w")
fout.write(fin.readline())
for line in fin:
    cvar = line.split("\t")
    cvar_id = cvar[0]
    cdna_var[cvar_id] = cvar
    if cvar_id not in cdna_vcf:
        hgvsp[cvar_id] = cvar
        fout.write(line)
    if cvar_id in long_var_id:
        f_umap_long_var.write(line)
fin.close()
fout.close()

f_umap_long_var.close()

fin = open(prot_var_fn)
fin.readline()
for line in fin:
    cvar = line.split("\t")
    cvar_id = cvar[0]
    hgvsp[cvar_id] = cvar
fin.close()


#Create mapping file for gDNA/cDNA civic info
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
fout_dna = open(os.path.join(out_dir, "civic_gdcmaf_mapping_dna.tsv"), "w")
fout_dna.write("civic_var_id\tcivic_gene_id\tsource\tchromosome\tstart_position" +
               "\treference_allele\talternative_allele\n")
fout_dna_full = open(os.path.join(out_dir, "civic_gdcmaf_mapping_dna_full_info.tsv"), "w")
fout_dna_full.write("civic_var_id\tcivic_gene_id\tsource\tchromosome\tstart_position" +
                    "\treference_allele\talternative_allele\ttranscript\ttranscript_2" +
                    "\tensembl_version\tref_build\tentrez_id\tentrez_name\tcivic_var_name" +
                    "\tcivic_var_types\tcivic_hgvs_exp\n")
for cvar_id, cvar in sorted(gdna_vcf.items(), key = lambda x:int(x[0])):
    fout_dna.write("\t".join(cvar) + "\n")
    v = gdna_var[cvar_id]
    fout_dna_full.write("\t".join(cvar + v[2:3] + v[4:7] + v[8:13]) + "\n")

for cvar_id, cvar in sorted(cdna_vcf.items(), key = lambda x:int(x[0])):
    fout_dna.write("\t".join(cvar) + "\n")
    v = cdna_var[cvar_id]
    fout_dna_full.write("\t".join(cvar + v[2:3] + v[4:7] + v[8:13]) + "\n")

fout_dna.close()
fout_dna_full.close()

#Create mapping file for prot civic info
fout_prot = open(os.path.join(out_dir, "civic_gdcmaf_mapping_prot.tsv"), "w")
fout_prot.write("civic_var_id\tcivic_gene_id\thugo_symbol\tgene\thgvs.p\tsource\n")

fout_prot_full = open(os.path.join(out_dir, "civic_gdcmaf_mapping_prot_full_info.tsv"), "w")
fout_prot_full.write("civic_var_id\tcivic_gene_id\thugo_symbol\tgene\thgvs.p\tsource" +
                     "\ttranscript\ttranscript_2\tensembl_version\tref_build\tentrez_id" +
                     "\tentrez_name\tcivic_var_name\tcivic_var_types\tcivic_hgvs_exp\n")
f_umap_g = open(os.path.join(out_dir, "parsed_civic_variants.unmapped_no_genecode.tsv"), "w")
f_umap_p = open(os.path.join(out_dir, "parsed_civic_variants.unmapped_no_p.tsv"), "w")
hp = hgvs.parser.Parser()
for civic_var_id, civic_var in sorted(hgvsp.items(), key = lambda x:int(x[0])):
    civic_gene_id = civic_var[7]
    civic_entrez_name = civic_var[9]
    gene = ""
    if civic_entrez_name in genecode:
        gene = genecode[civic_entrez_name]
    else:
        f_umap_g.write("\t".join(civic_var))
        continue
    aa1_p = []
    if civic_var[19] != "": #get hgvs p
        for p in civic_var[19].split(";"):
            var_p = hp.parse_hgvs_variant(p)
            var_p_aa1 = str(var_p.format(conf={"p_3_letter": False}))
            aa1_p.append(var_p_aa1.split(":")[1])
    elif civic_var[21] != "": #get vname p
        aa1_p.append("p." + civic_var[21])
    else:
        f_umap_p.write("\t".join(civic_var))
        continue
    for prot in aa1_p:
        fout_prot.write("\t".join([civic_var_id, civic_gene_id, civic_entrez_name, gene, prot, "Protein"]) + "\n")
        fout_prot_full.write("\t".join([civic_var_id, civic_gene_id, civic_entrez_name, gene, prot, "Protein"] \
                                       + civic_var[2:3] + civic_var[4:7] + civic_var[8:13]) + "\n")

fout_prot.close()
fout_prot_full.close()
