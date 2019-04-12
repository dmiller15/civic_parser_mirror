import os
import sys
import gzip

f1 = open(sys.argv[1])
f1.readline()
civic_var = {}
for line in f1:
    tmp = line.strip().split("\t")
    civic_var[tmp[7]] = tmp[0]
print(len(civic_var))

genecode = {}
for line in gzip.open(sys.argv[2]):
    if line[0] == "#":
        continue
    tmp = line.split("\t")
    annot = tmp[8].split(";")
    annot_dict = {}
    for tmp2 in annot:
    	tmp3 = tmp2.strip().split(" ")
    	if len(tmp3) >= 2:
	    	annot_dict[tmp3[0].strip()] = tmp3[1].strip('\"')
    try:
        genecode[annot_dict["gene_name"]] = annot_dict["gene_id"]
    except:
        print "Error", sys.exc_info()[0], line, annot_dict

print len(genecode)
'''
for k, v in civic_var.iteritems():
    if k not in genecode:
        print k, v
'''

cv = ["UGT1A", "MRE11"]
for k, v in genecode.iteritems():
    if cv[0] in k or k in cv[0]:
        print cv[0], civic_var[cv[0]], k, v
    if cv[1] in k or k in cv[1]:
        print cv[1], civic_var[cv[1]], k, v