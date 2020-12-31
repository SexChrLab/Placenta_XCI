# In this script, after running the script calc_median_allele_balance_per_tissue.py, we want to select only the tissues where there are at least 10 samples
import pandas as pd
import json

chrX_out_fn = open('tissues_with_more_than_10_samples_chrX.txt', 'w')
chr8_out_fn = open('tissues_with_more_than_10_samples_chr8.txt', 'w')

with open('/scratch/tphung3/Placenta_XCI/gtex/01_download_data/gtex_counts_analysis_config.json') as json_file:
    data = json.load(json_file)
    for tissue in sorted(data['tissues_ids']):
        tissue_data = pd. read_csv('/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/' + tissue + '_chrX_nonpars_median_allele_balance.txt', sep='\t')
        if tissue_data.shape[0] > 10:
            out = [tissue, str(tissue_data.shape[0])]
            print ('\t'.join(out), file=chrX_out_fn)
        else:
            print (tissue)

with open('/scratch/tphung3/Placenta_XCI/gtex/01_download_data/gtex_counts_analysis_config.json') as json_file:
    data = json.load(json_file)
    for tissue in sorted(data['tissues_ids']):
        tissue_data = pd. read_csv('/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/autosomes_WASP_raw_counts_allele_balance/chr8/' + tissue + '_chr8_median_allele_balance.txt', sep='\t')
        if tissue_data.shape[0] > 10:
            out = [tissue, str(tissue_data.shape[0])]
            print ('\t'.join(out), file=chr8_out_fn)
