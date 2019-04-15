import os
import sys
import hgvs.parser

gdc_variant_types_short = ["sub", "dup", "del", "ins", "delins", "fs"]
gdc_variant_types_long = ["point_mutation", "duplication", "small_deletion", "insertion", "delins", "frameshift"]

#To get genecode gene_name and gene_id
genecode = {}
f = open("data/gencode.gene.info.v22.tsv")
f.readline()
for line in f:
    tmp = line.split("\t")
    genecode[tmp[1]] = tmp[0]
print len(genecode)

#To get civic gene_id and var_id
civic_var = {}
f = open("data/civic_variants_extraction.tsv")
f.readline()
for line in f:
    tmp = line.strip().split("\t")
    civic_var[tmp[0]] = tmp[1]

#To get parsed hgvs.g to eliminate them from hgvs.p results
gDNA_var = {}
f = open("data/results/civic_gdcmaf_mapping_dna.tsv")
f.readline()
for line in f:
    tmp = line.strip().split("\t")
    gDNA_var[tmp[0]] = True


#Write result to file
fout = open(sys.argv[1], "w")
fout.write("civic_var_id\tcivic_gene_id\thugo_symbol\tgene\thgvs.p\tsource\n")

hp = hgvs.parser.Parser()
#To extract and process hgvs.p from all civic variants
fin = open("data/gdc_parsed_civic_variants_all.tsv")
fin.readline()
for line in fin:
    tmp = line.split("\t")
    civic_var_id = tmp[0]
    if civic_var_id in gDNA_var:
        continue
    civic_gene_id = civic_var[civic_var_id]
    civic_entrez_name = tmp[7]
    aa1_p = ""
    #Get vname p
    if tmp[20] != "" and tmp[21] in gdc_variant_types_long:
        aa1_p = "p." + tmp[19]
    #Get hgvs p
    elif tmp[17] != "" and tmp[18] in gdc_variant_types_short:
        var_p = hp.parse_hgvs_variant(tmp[17])
        var_p_aa1 = str(var_p.format(conf={"p_3_letter": False}))
        aa1_p = var_p_aa1.split(":")[1]
    else:
        print tmp[17:22]
        continue
    gene = ""
    if civic_entrez_name in genecode:
        gene = genecode[civic_entrez_name]
    else:
        print "non-geneconde-mapping", civic_entrez_name
    fout.write("\t".join([civic_var_id, civic_gene_id, civic_entrez_name, gene, aa1_p, "Protein"]) + "\n")

fout.close()
