# First, find the individuals where both of the sites are showing skewness
# Previously, 52/58 placenta samples are showing skewness and we use this for the gene escape analysis
# For showing heterogeneity within sample, we are considering just the samples with data for both sites
from collections import defaultdict
paired_samples = set()
with open('/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chrX/placenta_chrX_nonpars_median_allele_balance.txt', 'r') as f:
    samples = defaultdict(list)
    for line in f:
        if not line.startswith('tissue'):
            items = line.rstrip('\n').split('\t')
            if items[2] != "OBG0175":
                if float(items[7]) > 0.8:
                    id = items[2] + '_' + items[3]
                    samples[items[2]].append(id)

    for sample in samples:
        if len(samples[sample]) == 2:
            paired_samples.add(samples[sample][0])
            paired_samples.add(samples[sample][1])

# Second, find escape genes
escape_genes = set()
with open("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/chrX_escaping_samples_prop_per_gene.tsv", "r") as f:
    for line in f:
        if not line.startswith("gene"):
            items = line.rstrip("\n").split("\t")
            if items[5] == "Escape":
                escape_genes.add(items[0])

escape_outfile = open("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_s10/escape.csv", "w")

for i in escape_genes:
    with open("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/" + i + "_allele_balance_for_placenta_samples.tsv", "r") as f:
        for line in f:
            if not line.startswith("sampleID"):
                items = line.rstrip("\n").split("\t")
                if items[0] in paired_samples:
                    if items[1] != "NA":
                        if float(items[1]) > 0.8:
                            out = [items[0], "Inactivated", i]
                            print(",".join(out), file=escape_outfile)
                        else:
                            out = [items[0], "Escape", i]
                            print(",".join(out),
                                  file=escape_outfile)
                    else:
                        out = [items[0], "NA", i]
                        print(",".join(out),
                              file=escape_outfile)

# Third, find variable escape genes
variable_escape_genes = set()
with open("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/chrX_escaping_samples_prop_per_gene.tsv", "r") as f:
    for line in f:
        if not line.startswith("gene"):
            items = line.rstrip("\n").split("\t")
            if items[5] == "Variable_escape":
                variable_escape_genes.add(items[0])

variable_escape_outfile = open("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_s10/variable_escape.csv", "w")

for i in variable_escape_genes:
    with open("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/" + i + "_allele_balance_for_placenta_samples.tsv", "r") as f:
        for line in f:
            if not line.startswith("sampleID"):
                items = line.rstrip("\n").split("\t")
                if items[0] in paired_samples:
                    if items[1] != "NA":
                        if float(items[1]) > 0.8:
                            out = [items[0], "Inactivated", i]
                            print(",".join(out), file=variable_escape_outfile)
                        else:
                            out = [items[0], "Escape", i]
                            print(",".join(out),
                                  file=variable_escape_outfile)
                    else:
                        out = [items[0], "NA", i]
                        print(",".join(out),
                              file=variable_escape_outfile)

# Fourth, find inactivated genes
inactivated_genes = set()
with open("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/chrX_escaping_samples_prop_per_gene.tsv", "r") as f:
    for line in f:
        if not line.startswith("gene"):
            items = line.rstrip("\n").split("\t")
            if items[5] == "Inactivated":
                inactivated_genes.add(items[0])

inactivated_outfile = open("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_s10/inactivated.csv", "w")

for i in inactivated_genes:
    with open("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/" + i + "_allele_balance_for_placenta_samples.tsv", "r") as f:
        for line in f:
            if not line.startswith("sampleID"):
                items = line.rstrip("\n").split("\t")
                if items[0] in paired_samples:
                    if items[1] != "NA":
                        if float(items[1]) > 0.8:
                            out = [items[0], "Inactivated", i]
                            print(",".join(out), file=inactivated_outfile)
                        else:
                            out = [items[0], "Escape", i]
                            print(",".join(out),
                                  file=inactivated_outfile)
                    else:
                        out = [items[0], "NA", i]
                        print(",".join(out),
                              file=inactivated_outfile)