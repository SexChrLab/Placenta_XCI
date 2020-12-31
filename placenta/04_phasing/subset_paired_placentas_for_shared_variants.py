# This is the first script in phasing where we are subsetting the paired placentas for shared variants only
import pandas as pd
import json
import sys

chr = sys.argv[1]

totalCount_threshold = 10 #Change this if needs be

with open('/data/CEM/wilsonlab/projects/placenta/Placenta_XCI/placenta/03_analyze_ase/scripts/analyze_ase_config.json') as json_file:
    data = json.load(json_file)
    for sample_id in data['dna_rna_placenta']:
        # Read in the files
        site_A = pd.read_csv('/data/CEM/wilsonlab/projects/placenta/Placenta_XCI/placenta/03_analyze_ase/results/' + chr + '/' + sample_id + '_' + data['dna_rna_placenta'][sample_id][0] + '_' + chr + '_allele_balance.tsv', sep='\t')
        site_B = pd.read_csv('/data/CEM/wilsonlab/projects/placenta/Placenta_XCI/placenta/03_analyze_ase/results/' + chr + '/' + sample_id + '_' + data['dna_rna_placenta'][sample_id][1] + '_' + chr + '_allele_balance.tsv', sep='\t')

        # Find shared variants
        shared_variants = pd.merge(site_A, site_B, on='position', how='inner') #This returns the variants that are shared between site A and site B

        # Filter based on total_count threshold. This removes variants if the total count (for both site A and site B) is less than the threshold.
        shared_variants_threshold = shared_variants[(shared_variants['total_count_x']>totalCount_threshold) & (shared_variants['total_count_y']>totalCount_threshold)]

        shared_variants_threshold.to_csv('/data/CEM/wilsonlab/projects/placenta/Placenta_XCI/placenta/08_phasing/paired_placentas_shared_variants/' + sample_id + '_paired_placentas_shared_variants_' + chr + '.csv', index=False)

        # This print out: column 1: sample name, column 2: the number of variants in site A where allele balance is greater than 0.8, column 3: the number of variants in site B where allele balance is greater than 0.8.
        out = [sample_id, str(shared_variants_threshold.shape[0]), str(len(shared_variants_threshold[shared_variants_threshold['allele_balance_x']>0.8])), str(len(shared_variants_threshold[shared_variants_threshold['allele_balance_y']>0.8]))]
        print (','.join(out))
