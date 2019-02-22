import os
import sys
import re
from civicpy import civic
import hgvs.parser
import hgvs.dataproviders.uta
import hgvs.assemblymapper


hp = hgvs.parser.Parser()
hdp = hgvs.dataproviders.uta.connect()
am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh37', alt_aln_method='splign', replace_reference=True)

all_ids = range(4, 5)
all_var = civic.get_variants_by_ids(all_ids)

for v in all_var:
	for hgvs_var in v.hgvs_expressions:
		try:
			parsed_var = hp.parse_hgvs_variant(hgvs_var)
			t = parsed_var.posedit.edit.type
			if parsed_var.type == "g":
				print(t)
			elif parsed_var.type == "c":
				try:
					print("c2g parsed_var", am.c_to_g(parsed_var))
				except:
					print("unable to convert c to g for variant " + str(parsed_var))
			elif parsed_var.type == "p":
				print(t)
		except:
			print("unable to parse hgvs variant" + str(hgvs_var))
