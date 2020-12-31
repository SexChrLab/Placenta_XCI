# In this script, we want to compute: using the placenta data, for each gene, the proportion of samples that are escaping and that are inactivated
# Usage: python compute_escaping_samples_prop_per_gene.py chrX_escaping_samples_prop_per_gene.tsv

import sys
import os
import pandas as pd

threshold = 5

if os.path.exists(sys.argv[1]):
    outfile = open(sys.argv[1], "a")
else:
    outfile = open(sys.argv[1], "w")
    header = ["gene", "n_samples", "n_total", "prop", "status", "gene_status"]
    print ('\t'.join(header), file=outfile)

# Obtain a list of genes (690 genes)
with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/expressed_genes_all_samples.txt', 'r') as f:
    genes_list = [line.rstrip('\n') for line in f]


for gene in genes_list:
    placenta = pd.read_csv(
        '/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/' + gene + '_allele_balance_for_placenta_samples.tsv',
        sep='\t')
    placenta_noNA = placenta[placenta['allele_balance'].notnull()]
    placenta_noNA_count = placenta_noNA.shape[0]

    if placenta_noNA_count >= threshold:
        placenta_noNA_median = placenta_noNA['allele_balance'].median()
        placenta_noNA_inactivated = placenta_noNA[placenta_noNA['allele_balance'] >= 0.8]
        placenta_noNA_inactivated_count = placenta_noNA_inactivated.shape[0]
        placenta_noNA_inactivated_prop = placenta_noNA_inactivated_count / placenta_noNA_count

        out_1 = [gene, str(placenta_noNA_inactivated_count), str(placenta_noNA_count),
                 str(placenta_noNA_inactivated_prop), "Inactivated"]
        out_2 = [gene, str(placenta_noNA_count - placenta_noNA_inactivated_count), str(placenta_noNA_count),
                 str(1 - placenta_noNA_inactivated_prop), "Escape"]

        if placenta_noNA_inactivated_prop >= 0.7 and placenta_noNA_median >= 0.8:
            out_1.append('Inactivated')
            out_2.append('Inactivated')
        elif placenta_noNA_inactivated_prop <= 0.3 and placenta_noNA_median <=0.75:
            out_1.append('Escape')
            out_2.append('Escape')
        else:
            out_1.append('Variable_escape')
            out_2.append('Variable_escape')

        print("\t".join(out_1), file=outfile)
        print("\t".join(out_2), file=outfile)
