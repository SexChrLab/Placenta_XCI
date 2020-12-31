# In this script, I tabulate the number of individuals for each gtex tissue and the number of individuals per gtex tissue that is highly skewed

import json
from collections import defaultdict

out = defaultdict(list)

outfile = open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/tabulate_individuals.csv', 'w')
header = ['tissue', 'n_samples', 'n_samples_skewed']
print (','.join(header), file=outfile)

with open('/scratch/tphung3/Placenta_XCI/gtex/01_download_data/gtex_counts_analysis_config.json') as json_file:
    data = json.load(json_file)
    for tissue in data['tissues_ids']:
        out[tissue].append(len(data[tissue]))

with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/escape_genes_config.json') as json_file:
    data = json.load(json_file)
    for tissue in data['tissues_ids']:
        out[tissue].append(len(data[tissue]))

for k in sorted(out):
    if len(out[k]) == 1:
        to_print = [k, str(out[k][0]), str(0)]
        print (','.join(to_print), file=outfile)
    else:
        to_print = [k, str(out[k][0]), str(out[k][1])]
        print (','.join(to_print), file=outfile)
