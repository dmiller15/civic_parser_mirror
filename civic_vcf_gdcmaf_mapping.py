import os
import sys

gdc_variant_types_short = ["sub", "dup", "del", "ins", "delins", "fs"]
gdc_variant_types_long = ["point_mutation", "duplication", "small_deletion", "insertion", "delins", "frameshift"]

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

results_dir = sys.argv[1]
if not os.path.exists(os.path.join(results_dir, "mapping")):
    os.makedirs(os.path.join(results_dir, "mapping"))

fout = open(os.path.join(results_dir, "mapping", "civic_gdcmaf_mapping_dna.tsv", "w"))
fout.write("civic_var_id\tcivic_gene_id\tsource\tchromosome\tstart_position\treference_allele\talternative_allele\n")

fin = open(os.path.join(results_dir, "gDNA", "gDNA_transvar_parsed_civic_variants_combined_lifted_over.vcf"))
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

fin = open(os.path.join(results_dir, "cDNA", "cDNA_transvar_parsed_civic_variants_combined_lifted_over.vcf"))
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
