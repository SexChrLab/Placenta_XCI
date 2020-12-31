# In this script, we are adding the genes to the config file (escape_genes_config.json)
# Previously we used the script `find_samples_skewed.py` to generate the config file escape_genes_config.json
# Now we are adding the genes

import json

with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/escape_genes_config.json') as json_file:
    data = json.load(json_file)

data['genes'] = []
with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/expressed_genes_all_samples.txt', 'r') as f:
    for line in f:
        data['genes'].append(line.rstrip('\n'))

with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/escape_genes_config.json', 'w') as f:
    json.dump(data, f)
