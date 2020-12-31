# In this script, we are generating for each gene
# A file for GTEX
# A file for placenta
# Comment out deciduas

import json

with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/escape_genes_config.json') as json_file:
    data = json.load(json_file)
    gtex_samples = data["gtex_samples"]
    placenta_samples = data["full_term_placenta"]
    # decidua_females_samples = data["decidua_females"]
    # decidua_males_samples = data["decidua_males"]

gtex_dict = {}
full_term_placenta_dict = {}
# decidua_females_dict = {}
# decidua_males_dict = {}

with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/expressed_genes_all_samples.txt', 'r') as f:
    for line in f:
        gtex_dict[line.rstrip('\n')] = {}
        full_term_placenta_dict[line.rstrip('\n')] = {}
        # decidua_females_dict[line.rstrip('\n')] = {}
        # decidua_males_dict[line.rstrip('\n')] = {}
        for sample in gtex_samples:
            gtex_dict[line.rstrip('\n')][sample] = 'NA'
        for sample in placenta_samples:
            full_term_placenta_dict[line.rstrip('\n')][sample] = 'NA'
        # for sample in decidua_females_samples:
        #     decidua_females_dict[line.rstrip('\n')][sample] = 'NA'
        # for sample in decidua_males_samples:
        #     decidua_males_dict[line.rstrip('\n')][sample] = 'NA'

for sample in gtex_samples:
    with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/asereadcounter_geneinfo/chrX/' + sample + '_chrX_allele_balance_per_gene.tsv', 'r') as f:
        for line in f:
            items = line.rstrip('\n').split('\t')
            gtex_dict[items[0]][sample] = items[1]

for sample in placenta_samples:
    with open('/scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/asereadcounter_geneinfo/chrX/' + sample + '_chrX_allele_balance_per_gene.tsv', 'r') as f:
        for line in f:
            items = line.rstrip('\n').split('\t')
            full_term_placenta_dict[items[0]][sample] = items[1]

# for sample in decidua_females_samples:
#     with open('/scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/asereadcounter_geneinfo/chrX/' + sample + '_chrX_allele_balance_per_gene.tsv', 'r') as f:
#         for line in f:
#             items = line.rstrip('\n').split('\t')
#             decidua_females_dict[items[0]][sample] = items[1]
#
# for sample in decidua_males_samples:
#     with open('/scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/asereadcounter_geneinfo/chrX/' + sample + '_chrX_allele_balance_per_gene.tsv', 'r') as f:
#         for line in f:
#             items = line.rstrip('\n').split('\t')
#             decidua_males_dict[items[0]][sample] = items[1]

for gene in gtex_dict:
    outfile = open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/' + gene + '_allele_balance_for_gtex_samples.tsv', 'w')
    header = ['sampleID', 'allele_balance']
    print ('\t'.join(header), file=outfile)
    for sample in gtex_dict[gene]:
        out = [sample, gtex_dict[gene][sample]]
        print ('\t'.join(out), file=outfile)

for gene in full_term_placenta_dict:
    outfile = open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/' + gene + '_allele_balance_for_placenta_samples.tsv', 'w')
    header = ['sampleID', 'allele_balance']
    print ('\t'.join(header), file=outfile)
    for sample in full_term_placenta_dict[gene]:
        out = [sample, full_term_placenta_dict[gene][sample]]
        print ('\t'.join(out), file=outfile)

# for gene in decidua_females_dict:
#     outfile = open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/' + gene + '_allele_balance_for_decidua_females_samples.tsv', 'w')
#     header = ['sampleID', 'allele_balance']
#     print ('\t'.join(header), file=outfile)
#     for sample in decidua_females_dict[gene]:
#         out = [sample, decidua_females_dict[gene][sample]]
#         print ('\t'.join(out), file=outfile)
#
# for gene in decidua_males_dict:
#     outfile = open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/' + gene + '_allele_balance_for_decidua_males_samples.tsv', 'w')
#     header = ['sampleID', 'allele_balance']
#     print ('\t'.join(header), file=outfile)
#     for sample in decidua_males_dict[gene]:
#         out = [sample, decidua_males_dict[gene][sample]]
#         print ('\t'.join(out), file=outfile)
