import os
import sys
import gzip

civic_transvar = {}
for line in open(sys.argv[1]):
    if line[0] == "#":
        continue
    tmp = line.strip().split("\t")
    civic_transvar[tmp[0] + ":"  + tmp[1] + ":" + tmp[3] + ":" + tmp[4]] = 0
print(len(civic_transvar))

for filename in os.listdir(sys.argv[2]):
    tmp = filename.split(".")
    if tmp[0] == "TCGA" and tmp[2] == "mutect" and tmp[6] == "protected":
        print(os.path.join(sys.argv[2], filename))
        for line in gzip.open(os.path.join(sys.argv[2], filename)):
            if line[0] == "#" or "Hugo_Symbol" in line:
                continue
            tmp = line.strip().split("\t")
            vcf_region = tmp[121].strip().split(":")
            if len(vcf_region) == 5:
                k = vcf_region[0] + ":"  + vcf_region[1] + ":" + vcf_region[3] + ":" + vcf_region[4]
                if k in civic_transvar:
                    civic_transvar[k] += 1
            else:
                print("abnormal vcf_region", vcf_region)

f = open(sys.argv[3], "w")
for k, v in civic_transvar.items():
    f.write(str(k) + "\t" + str(v) + "\n")
f.close()
