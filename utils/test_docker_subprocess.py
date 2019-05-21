import sys
from subprocess import check_output

try:
    parsed_var = "ENST00000389048.3:c.3824G>A"
    print("using TransVar", str(parsed_var))
    #out = check_output(["ls", "-l"])
    print("TransVar cmd", " ".join(["docker", "run", "-v", "/mnt/tools/transvar.download/:/ref", "-t", \
        "quay.io/ncigdc/transvar:2.5.3.hg19", "transvar", "canno", "--ensembl", "--reference", "/ref/hg19.fa", \
        "--gseq", "-i", "\"" + str(parsed_var) + "\""]))
    out = check_output(["docker", "run", "-v", "/mnt/tools/transvar.download/:/ref", "-t", \
                        "quay.io/ncigdc/transvar:2.5.3.hg19", "transvar", "canno", "--ensembl", "--reference", "/ref/hg19.fa", \
                        "--gseq", "-i", shlex.split("\"" + str(parsed_var) + "\"")])
    a = out.decode('utf-8').strip()
    for l in a.split("\n"):
        print(l)

except:
    print("TransVar c2g conversion error", sys.exc_info()[0])
    
