# In this script, we are returning genes that have at least one heterozygous and expressed variant that are found across all samples (after threshold filtering)
# This is merging genes from all samples

import os
import json

gtex_files_dir = '/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/asereadcounter_geneinfo/chrX/'
placenta_files_dir = '/scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/asereadcounter_geneinfo/chrX/' #this is because I processed this previously in wes_genotyping
threshold = 10

all_genes = set()

# gtex
for file in [f for f in os.listdir(gtex_files_dir) if f.endswith('allele_balance_geneinfo_unique.tsv')]:
    with open(os.path.join(gtex_files_dir, file), 'r') as f:
        for line in f:
            items = line.rstrip('\n').split('\t')
            if float(items[11]) > threshold:
                all_genes.add(items[3])

# placenta
placenta_ids = set()
with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/escape_genes_config.json') as json_file:
    data = json.load(json_file)
    for i in data["full_term_placenta"]:
        placenta_ids.add(i.split("_")[0])

for file in os.listdir(placenta_files_dir):
    if file.endswith('allele_balance_geneinfo_unique.tsv'):
        if file.split("_")[0] in placenta_ids:
            with open(os.path.join(placenta_files_dir, file), 'r') as f:
                for line in f:
                    items = line.rstrip('\n').split('\t')
                    if float(items[11]) > threshold:
                        all_genes.add(items[3])

outfile = open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/expressed_genes_all_samples.txt', 'w')
for i in all_genes:
    print (i, file=outfile)