# In this script, I want to subset the ASEReadCounter counts downloaded from GTEX for each rna_id
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Subset gtex counts for each rna_id')
parser.add_argument('--rna_id',required=True,help='Input the rna_id')
parser.add_argument('--counts_dir',required=True,help='Input the directory to the counts. For examples: /data/mwilsons/public_data/controlled_access/gtex/version8/chrX_WASP_raw_counts/')
parser.add_argument('--ending_str',required=True,help='For examples: .v8.readcounts.chrX.txt.gz')
parser.add_argument('--chr',required=True,help='Input either chrX or autosomes')
parser.add_argument('--outfile',required=True,help='Input the path to the outfile')

args = parser.parse_args()

counts_dir = args.counts_dir

rna_id = args.rna_id
j = rna_id.split('-')
sample_id = j[0] + '-' + j[1]
data = pd.read_csv(os.path.join(counts_dir, sample_id + args.ending_str), sep='\t', compression='gzip')
if args.chr == 'chrX':
    data_subset = data[data['sampleID'] == rna_id]
    data_subset.to_csv(args.outfile, sep='\t', index=False)
else:
    data_subset = data[data['SAMPLE_ID'] == rna_id]
    data_subset.to_csv(args.outfile, sep='\t', index=False)
