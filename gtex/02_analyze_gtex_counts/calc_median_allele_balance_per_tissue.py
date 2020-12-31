import json
import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Calculate median allele balance per tissue.')
parser.add_argument('--tissue',required=True,help='Input the name of the tissue')
parser.add_argument('--chromosome',required=True,help='Input the chromosome number')
parser.add_argument('--files_dir',required=True,help='Input the path to the directory where the tsv file of the allele balance is. For example: /scratch/tphung3/PlacentaSexDiff/F_gtex_counts/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/')
parser.add_argument('--threshold',required=True,help='Input the threshold to subset the total count with')
parser.add_argument('--include_pars',required=True,help='Yes if PARs are included in the calculation. No otherwise.')
parser.add_argument('--out',required=True,help='Input the path to the output file')

args = parser.parse_args()

tissue = args.tissue
chromosome = args.chromosome
threshold = float(args.threshold)
files_dir = args.files_dir
out = args.out

if not os.path.exists(out):
    header = ['tissue', 'chr', 'rna_id', 'num_var', 'median_allele_balance', 'num_var_subset', 'median_allele_balance_subset']
    outfile = open(out, 'w')
    print ('\t'.join(header), file=outfile)
else:
    outfile = open(out, 'a')

with open('/scratch/tphung3/PlacentaSexDiff/F_gtex_counts/01_download_data/gtex_counts_analysis_config.json') as json_file:
    data = json.load(json_file)
    for rna_id in data[tissue]:
        data = pd.read_csv(os.path.join(files_dir, rna_id + '_allele_balance.tsv'), sep='\t') #header of input file is: header = ['chr', 'position', 'ref_allele', 'alt_allele', 'ref_count', 'alt_count', 'total_count', 'allele_balance']
        if args.include_pars == 'yes':
            data_nrow = len(data.index)
            median_allele_balance = data['allele_balance'].median()

            # Subset based on threshold of total count
            data_subset = data[data['total_count'] > threshold]
            data_subset_nrow = len(data_subset.index)
            median_allele_balance_subset = data_subset['allele_balance'].median()

            output = [tissue, chromosome, rna_id, str(data_nrow), str(median_allele_balance), str(data_subset_nrow), str(median_allele_balance_subset)]
            print ('\t'.join(output), file=outfile)
        else: #remove pars
            data_rmpars_1 = data[(data['position'] < 10001) | (data['position'] > 2781479)]
            data_rmpars_both = data_rmpars_1[(data_rmpars_1['position'] < 155701383) | (data_rmpars_1['position'] > 156030895)]
            data_rmpars_both_nrow = len(data_rmpars_both.index)
            median_allele_balance = data_rmpars_both['allele_balance'].median()

            # Subset based on threshold of total count
            data_rmpars_both_subset = data_rmpars_both[data_rmpars_both['total_count'] > threshold]
            data_rmpars_both_subset_nrow = len(data_rmpars_both_subset.index)
            median_allele_balance_subset = data_rmpars_both_subset['allele_balance'].median()

            output = [tissue, chromosome, rna_id, str(data_rmpars_both_nrow), str(median_allele_balance), str(data_rmpars_both_subset_nrow), str(median_allele_balance_subset)]
            print ('\t'.join(output), file=outfile)

