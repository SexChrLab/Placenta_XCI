# In this script we are computing allele balance using the phased data
import random
import numpy as np
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(description='Compute allele balance using the phased data.')
parser.add_argument('--sample_id',required=True,help='Input the sample id. For example: OBG0044')
parser.add_argument('--chromosome',required=True,help='Input the chromosome. For example: chr8 or chrX')
parser.add_argument('--shared_variants',required=True,help='Input the path to the shared variants file. For example: /scratch/tphung3/PlacentaSexDiff/A_placenta/08_phasing/paired_placentas_shared_variants/OBG0044_paired_placentas_shared_variants_chrX.csv')
parser.add_argument('--out_summary',required=True,help='Input the outfile summarizing allele balance')
parser.add_argument('--out_phased_data',required=True,help='Input the outfile detailing the phased information')

args = parser.parse_args()

if os.path.exists(args.out_summary):
    out_summary = open(args.out_summary, 'a')
else:
    out_summary = open(args.out_summary, 'w')
    header = ['sample_id', 'site_A_phased_allele_balance', 'site_B_phased_allele_balance']
    print ('\t'.join(header), file=out_summary)

out_phased_data = open(args.out_phased_data, 'w')
header = ['chr_x', 'position', 'ref_allele_x', 'alt_allele_x', 'ref_count_x', 'alt_count_x', 'total_count_x', 'allele_balance_x', 'chr_y', 'ref_allele_y,alt_allele_y', 'ref_count_y', 'alt_count_y', 'total_count_y', 'allele_balance_y', 'phased_x', 'phased_y']
print (','.join(header), file=out_phased_data)

shared_variants_threshold = pd.read_csv(args.shared_variants)

site_A_phased_allele_balance_all = []
site_B_phased_allele_balance_all = []

a = len(shared_variants_threshold[shared_variants_threshold['allele_balance_x']>0.8])
b = len(shared_variants_threshold[shared_variants_threshold['allele_balance_y']>0.8])

if a > b:
    with open(args.shared_variants, 'r') as f:
        for line in f:
            if line.startswith(args.chromosome):
                items = line.rstrip('\n').split(',')
                ref_ratio = float(items[4])/float(items[6])
                alt_ratio = float(items[5])/float(items[6])
                if ref_ratio > 0.8:
                    site_A_phased_allele_balance = ref_ratio
                    site_B_phased_allele_balance = float(items[11])/float(items[13])
                    site_A_phased_allele_balance_all.append(site_A_phased_allele_balance)
                    site_B_phased_allele_balance_all.append(site_B_phased_allele_balance)

                    # for saving to out_phased_data
                    items.append(str(site_A_phased_allele_balance))
                    items.append(str(site_B_phased_allele_balance))
                    print(','.join(items), file=out_phased_data)

                elif alt_ratio > 0.8:
                    site_A_phased_allele_balance = alt_ratio
                    site_B_phased_allele_balance = float(items[12])/float(items[13])
                    site_A_phased_allele_balance_all.append(site_A_phased_allele_balance)
                    site_B_phased_allele_balance_all.append(site_B_phased_allele_balance)

                    # for saving to out_phased_data
                    items.append(str(site_A_phased_allele_balance))
                    items.append(str(site_B_phased_allele_balance))
                    print(','.join(items), file=out_phased_data)

                else:
                    allele_to_phase = random.choice(['ref', 'alt'])
                    if allele_to_phase == 'ref':
                        site_A_phased_allele_balance = ref_ratio
                        site_B_phased_allele_balance = float(items[11]) / float(items[13])
                        site_A_phased_allele_balance_all.append(site_A_phased_allele_balance)
                        site_B_phased_allele_balance_all.append(site_B_phased_allele_balance)

                        # for saving to out_phased_data
                        items.append(str(site_A_phased_allele_balance))
                        items.append(str(site_B_phased_allele_balance))
                        print(','.join(items), file=out_phased_data)

                    else:
                        site_A_phased_allele_balance = alt_ratio
                        site_B_phased_allele_balance = float(items[12]) / float(items[13])
                        site_A_phased_allele_balance_all.append(site_A_phased_allele_balance)
                        site_B_phased_allele_balance_all.append(site_B_phased_allele_balance)

                        # for saving to out_phased_data
                        items.append(str(site_A_phased_allele_balance))
                        items.append(str(site_B_phased_allele_balance))
                        print(','.join(items), file=out_phased_data)

else:
    with open(args.shared_variants, 'r') as f:
        for line in f:
            if line.startswith(args.chromosome):
                items = line.rstrip('\n').split(',')
                ref_ratio = float(items[11]) / float(items[13])
                alt_ratio = float(items[12]) / float(items[13])
                if ref_ratio > 0.8:
                    site_B_phased_allele_balance = ref_ratio
                    site_A_phased_allele_balance = float(items[4]) / float(items[6])
                    site_A_phased_allele_balance_all.append(site_A_phased_allele_balance)
                    site_B_phased_allele_balance_all.append(site_B_phased_allele_balance)

                    # for saving to out_phased_data
                    items.append(str(site_A_phased_allele_balance))
                    items.append(str(site_B_phased_allele_balance))
                    print(','.join(items), file=out_phased_data)

                elif alt_ratio > 0.8:
                    site_B_phased_allele_balance = alt_ratio
                    site_A_phased_allele_balance = float(items[5]) / float(items[6])
                    site_A_phased_allele_balance_all.append(site_A_phased_allele_balance)
                    site_B_phased_allele_balance_all.append(site_B_phased_allele_balance)

                    # for saving to out_phased_data
                    items.append(str(site_A_phased_allele_balance))
                    items.append(str(site_B_phased_allele_balance))
                    print(','.join(items), file=out_phased_data)

                else:
                    allele_to_phase = random.choice(['ref', 'alt'])
                    if allele_to_phase == 'ref':
                        site_B_phased_allele_balance = ref_ratio
                        site_A_phased_allele_balance = float(items[4]) / float(items[6])
                        site_A_phased_allele_balance_all.append(site_A_phased_allele_balance)
                        site_B_phased_allele_balance_all.append(site_B_phased_allele_balance)

                        # for saving to out_phased_data
                        items.append(str(site_A_phased_allele_balance))
                        items.append(str(site_B_phased_allele_balance))
                        print(','.join(items), file=out_phased_data)

                    else:
                        site_B_phased_allele_balance = alt_ratio
                        site_A_phased_allele_balance = float(items[5]) / float(items[6])
                        site_A_phased_allele_balance_all.append(site_A_phased_allele_balance)
                        site_B_phased_allele_balance_all.append(site_B_phased_allele_balance)

                        # for saving to out_phased_data
                        items.append(str(site_A_phased_allele_balance))
                        items.append(str(site_B_phased_allele_balance))
                        print(','.join(items), file=out_phased_data)

if np.median(site_A_phased_allele_balance_all) > np.median(site_B_phased_allele_balance_all):
    out = [args.sample_id, str(np.median(site_A_phased_allele_balance_all)), str(np.median(site_B_phased_allele_balance_all))]
    print ('\t'.join(out), file=out_summary)
else:
    out = [args.sample_id, str(np.median(site_B_phased_allele_balance_all)), str(np.median(site_A_phased_allele_balance_all))]
    print ('\t'.join(out), file=out_summary)