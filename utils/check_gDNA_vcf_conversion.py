import os
import sys

vcf_ori = {}
for line in open("data/gDNA/gDNA_parsed_civic_variants_combined.vcf"):
	tmp = line.split("\t")
	vcf_ori[":".join(tmp[0:5])] = tmp[5:]

vcf_fix = {}
for line in open("data/gDNA/gDNA_parsed_civic_variants_combined_fix.vcf"):
	tmp = line.split("\t")
	vcf_fix[":".join(tmp[0:5])] = tmp[5:]

for k, v in vcf_ori.iteritems():
	if k not in vcf_fix:
		print "ori not in fix", k, v
	else:
		if v != vcf_fix[k]:
			print "ori in fix, not equal", k, v, vcf_fix[k]
