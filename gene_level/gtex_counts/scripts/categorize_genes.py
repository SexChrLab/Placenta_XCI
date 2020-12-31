# In this script, I want to categorize genes based on the following threshold:
# 1. For a gene to be considered, it has to have at least 5 samples
# 2. Inactivated: >=70% of the samples where median allele balance >= 0.8
# 3. Escape: <=30% of the samples where median allele balance >=0.8
# 4. Variable escape: else
import pandas as pd

threshold = 5

# Obtain a list of genes (689 genes)
with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/expressed_genes_all_samples.txt', 'r') as f:
    genes_list = [line.rstrip('\n') for line in f]

# Set up the output file
outfile = open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/categorize_genes_gtex_placenta.csv', 'w')
header = ['Gene', 'GTEx_status', 'Placenta_status']
print (','.join(header), file=outfile)

for gene in sorted(genes_list):
    out = [gene]
    # gtex
    gtex = pd.read_csv('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/' + gene + '_allele_balance_for_gtex_samples.tsv', sep='\t')
    gtex_noNA = gtex[gtex['allele_balance'].notnull()]
    gtex_noNA_count = gtex_noNA.shape[0]

    if gtex_noNA_count >= threshold:
        gtex_noNA_median = gtex_noNA['allele_balance'].median()
        gtex_noNA_inactivated = gtex_noNA[gtex_noNA['allele_balance'] >= 0.8]
        gtex_noNA_inactivated_count = gtex_noNA_inactivated.shape[0]
        gtex_noNA_inactivated_prop = gtex_noNA_inactivated_count/gtex_noNA_count
        if gtex_noNA_inactivated_prop >= 0.7 and gtex_noNA_median >=0.8:
            out.append('Inactivated')
        elif gtex_noNA_inactivated_prop <= 0.3 and gtex_noNA_median <=0.75:
            out.append('Escape')
        else:
            out.append('Variable_escape')
    else:
        out.append('NA')

    # placenta
    placenta = pd.read_csv('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/' + gene + '_allele_balance_for_placenta_samples.tsv', sep='\t')
    placenta_noNA = placenta[placenta['allele_balance'].notnull()]
    placenta_noNA_count = placenta_noNA.shape[0]

    if placenta_noNA_count >= threshold:
        placenta_noNA_median = placenta_noNA['allele_balance'].median()
        placenta_noNA_inactivated = placenta_noNA[placenta_noNA['allele_balance'] >= 0.8]
        placenta_noNA_inactivated_count = placenta_noNA_inactivated.shape[0]
        placenta_noNA_inactivated_prop = placenta_noNA_inactivated_count/placenta_noNA_count
        if placenta_noNA_inactivated_prop >= 0.7 and placenta_noNA_median >= 0.8:
            out.append('Inactivated')
        elif placenta_noNA_inactivated_prop <= 0.3 and placenta_noNA_median <=0.75:
            out.append('Escape')
        else:
            out.append('Variable_escape')
    else:
        out.append('NA')

    # # Decidua females
    # decidua_females = pd.read_csv('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/' + gene + '_allele_balance_for_decidua_females_samples.tsv', sep='\t')
    # decidua_females_noNA = decidua_females[decidua_females['allele_balance'].notnull()]
    # decidua_females_noNA_count = decidua_females_noNA.shape[0]
    # if decidua_females_noNA_count >= threshold:
    #     decidua_females_noNA_median = decidua_females_noNA['allele_balance'].median()
    #     decidua_females_noNA_inactivated = decidua_females_noNA[decidua_females_noNA['allele_balance'] >= 0.8]
    #     decidua_females_noNA_inactivated_count = decidua_females_noNA_inactivated.shape[0]
    #     decidua_females_noNA_inactivated_prop = decidua_females_noNA_inactivated_count/decidua_females_noNA_count
    #     if decidua_females_noNA_inactivated_prop >= 0.7 and decidua_females_noNA_median >= 0.8:
    #         out.append('Inactivated')
    #     elif decidua_females_noNA_inactivated_prop <= 0.3 and decidua_females_noNA_median <= 0.75:
    #         out.append('Escape')
    #     else:
    #         out.append('Variable_escape')
    # else:
    #     out.append('NA')
    #
    # # Decidua males
    # decidua_males = pd.read_csv('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/' + gene + '_allele_balance_for_decidua_males_samples.tsv', sep='\t')
    # decidua_males_noNA = decidua_males[decidua_males['allele_balance'].notnull()]
    # decidua_males_noNA_count = decidua_males_noNA.shape[0]
    # if decidua_males_noNA_count >= threshold:
    #     decidua_males_noNA_median = decidua_males_noNA['allele_balance'].median()
    #     decidua_males_noNA_inactivated = decidua_males_noNA[decidua_males_noNA['allele_balance'] >= 0.8]
    #     decidua_males_noNA_inactivated_count = decidua_males_noNA_inactivated.shape[0]
    #     decidua_males_noNA_inactivated_prop = decidua_males_noNA_inactivated_count/decidua_males_noNA_count
    #     if decidua_males_noNA_inactivated_prop >= 0.7 and decidua_males_noNA_median >= 0.8:
    #         out.append('Inactivated')
    #     elif decidua_males_noNA_inactivated_prop <= 0.3 and decidua_males_noNA_median <= 0.75:
    #         out.append('Escape')
    #     else:
    #         out.append('Variable_escape')
    # else:
    #     out.append('NA')
    if out[1] != "NA" and out[2] != "NA":
        print (','.join(out), file=outfile)

# Sort for plotting heatmap
outfile_sorted = open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/categorize_genes_gtex_placenta_sorted.csv', 'w')

inactivated_both = []
escape_both = []
variable_both = []
gtex_inactivated_placenta_escape = []
gtex_escape_placenta_inactivated = []
gtex_inactivated_placenta_variable = []
gtex_escape_placenta_variable = []
gtex_variable_placenta_inactivated = []
gtex_variable_placenta_escape = []

with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/categorize_genes_gtex_placenta.csv', 'r') as f:
    for line in f:
        if line.startswith('Gene'):
            print (line.rstrip('\n'), file=outfile_sorted)
        else:
            items = line.rstrip('\n').split(',')
            if items[1] == 'Inactivated' and items[2] == 'Inactivated':
                inactivated_both.append(line.rstrip('\n'))
            elif items[1] == 'Escape' and items[2] == 'Escape':
                escape_both.append(line.rstrip('\n'))
            elif items[1] == 'Variable_escape' and items[2] == 'Variable_escape':
                variable_both.append(line.rstrip('\n'))
            elif items[1] == 'Inactivated' and items[2] == 'Escape':
                gtex_inactivated_placenta_escape.append(line.rstrip('\n'))
            elif items[1] == 'Escape' and items[2] == 'Inactivated':
                gtex_escape_placenta_inactivated.append(line.rstrip('\n'))
            elif items[1] == 'Inactivated' and items[2] == 'Variable_escape':
                gtex_inactivated_placenta_variable.append(line.rstrip('\n'))
            elif items[1] == 'Escape' and items[2] == 'Variable_escape':
                gtex_escape_placenta_variable.append(line.rstrip('\n'))
            elif items[1] == 'Variable_escape' and items[2] == 'Inactivated':
                gtex_variable_placenta_inactivated.append(line.rstrip('\n'))
            elif items[1] == 'Variable_escape' and items[2] == 'Escape':
                gtex_variable_placenta_escape.append(line.rstrip('\n'))


for i in inactivated_both:
    print (i, file=outfile_sorted)
for i in escape_both:
    print (i, file=outfile_sorted)
for i in variable_both:
    print (i, file=outfile_sorted)
for i in gtex_inactivated_placenta_escape:
    print (i, file=outfile_sorted)
for i in gtex_escape_placenta_inactivated:
    print (i, file=outfile_sorted)
for i in gtex_inactivated_placenta_variable:
    print (i, file=outfile_sorted)
for i in gtex_variable_placenta_inactivated:
    print (i, file=outfile_sorted)
for i in gtex_escape_placenta_variable:
    print (i, file=outfile_sorted)
for i in gtex_variable_placenta_escape:
    print (i, file=outfile_sorted)