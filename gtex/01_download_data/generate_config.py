# In this script, I will generate a config json file for use for these analyses
import json
from collections import defaultdict

# Setting up
data = {}
data['all_rna_ids'] = []
data['tissues_ids'] = []
tissues_ids = defaultdict(list)

tissues_ids_set = set()

# Remove the files that failed download
failed_files = set()
with open('/scratch/tphung3/Placenta_XCI/gtex/01_download_data/failed_files.txt', 'r') as f:
    for line in f:
        failed_files.add(line.rstrip('\n'))

# Remove males
females = set()
with open('/scratch/tphung3/Placenta_XCI/gtex/01_download_data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', 'r') as f:
    for line in f:
        if not line.startswith('SUBJID'):
            items = line.rstrip('\n').split('\t')
            if items[1] == '2':
                females.add(items[0])

# Processing
with open('/scratch/tphung3/Placenta_XCI/gtex/01_download_data/sample.tsv', 'r') as f:
    for line in f:
        if not line.startswith('entity'):
            items = line.rstrip('\n').split('\t')
            j = items[0].split('-')
            id = j[0] + '-' + j[1]
            if id not in failed_files:
                if id in females:
                    data['all_rna_ids'].append(items[0])
                    tissues_ids_set.add(items[4])
                    tissues_ids[items[4]].append(items[0])

for tissue_id in tissues_ids_set:
    data['tissues_ids'].append(tissue_id)

for tissue_id in tissues_ids:
    data[tissue_id] = tissues_ids[tissue_id]

with open('/scratch/tphung3/Placenta_XCI/gtex/01_download_data/gtex_counts_analysis_config.json', 'w') as outfile:
    json.dump(data, outfile)

print (len(data['all_rna_ids']))
print (len(tissues_ids_set))