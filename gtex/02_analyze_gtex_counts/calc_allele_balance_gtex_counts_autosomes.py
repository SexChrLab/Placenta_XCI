# This script calculates the allele balance
import argparse

parser = argparse.ArgumentParser(description='Calculate allele balance (unphased) using gtex counts for the autosomes.\n'
                                             'This is because the format of gtex counts for the autosomes is different from that for chrX.')
parser.add_argument('--input',required=True,help='Input the path to the output from asereadcounter')
parser.add_argument('--chr',required=True,help='Input chromosome number. For example: chr8')
parser.add_argument('--output',required=True,help='Input the path to the output file')

args = parser.parse_args()

# Header of input file:
# CHR     POS     VARIANT_ID      REF_ALLELE      ALT_ALLELE      SAMPLE_ID       SUBJECT_ID      TISSUE_ID       REF_COUNT       ALT_COUNT       TOTAL_COUNT     REF_RATIO       OTHER_ALLELE_COUNT      NULL_RATIO      BINOM_P BINOM_P_ADJUSTED        MAMBA_POST_SINGLETIS    MAMBA_POST_MULTITIS     GENOTYPE        VARIANT_ANNOTATION      GENE_ID LOW_MAPABILITY  MAPPING_BIAS_SIM
#         GENOTYPE_WARNING
input = args.input

# Format of output file:
output = open(args.output, 'w')
header = ['chr', 'position', 'ref_allele', 'alt_allele', 'ref_count', 'alt_count', 'total_count', 'allele_balance']
print ('\t'.join(header), file=output)

with open(input, 'r') as f:
    for line in f:
        if line.startswith(args.chr):
            items = line.rstrip('\n').split('\t')
            out = [items[0], items[1], items[3], items[4], items[8], items[9], items[10]]
            ref_ratio = float(items[8])/float(items[10])
            alt_ratio = float(items[9])/float(items[10])

            if ref_ratio > alt_ratio:
                out.append(str(ref_ratio))
            else:
                out.append(str(alt_ratio))
            print ('\t'.join(out), file=output)

