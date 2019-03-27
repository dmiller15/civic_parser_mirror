import os
import sys
import json
import requests

url = 'https://civicdb.org/api/variants?count=99999'
response = requests.get(url)
json_data = json.loads(response.text)

f = open(sys.argv[1], "w")
f.write("civic_var_id\tgen_id\tentrez_id\tentrez_name\n")
for x in json_data['records']:
    f.write(str(x['id']) + "\t" + str(x['gene_id']) + "\t" + str(x['entrez_id']) + "\t" + str(x['entrez_name']) + '\n')
f.close()