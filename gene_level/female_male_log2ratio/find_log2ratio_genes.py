# In this script, we want to return the female to male log2ratio for each category of genes (escape, inactivated, and variable escape)
log2ratio = {}
with open("/scratch/tphung3/Placenta_XCI/gene_level/female_male_log2ratio/DEGs_placentas_genesFvsM_fdr1_lfc0_chrX.txt", "r") as f:
    for line in f:
        if not line.startswith("row"):
            cols = line.rstrip("\n").split("\t")
            log2ratio[cols[1]] = cols[4]

inactivated_out = open("/scratch/tphung3/Placenta_XCI/gene_level/female_male_log2ratio/inactivated_genes_log2ratio.tsv", "w")
escape_out = open("/scratch/tphung3/Placenta_XCI/gene_level/female_male_log2ratio/escape_genes_log2ratio.tsv", "w")
variable_out = open("/scratch/tphung3/Placenta_XCI/gene_level/female_male_log2ratio/variable_genes_log2ratio.tsv", "w")

with open("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/categorize_genes_gtex_placenta_sorted.csv", "r") as f:
    for line in f:
        if not line.startswith("Gene"):
            cols = line.rstrip("\n").split(",")
            if cols[2] == "Inactivated":
                if cols[0] in log2ratio:
                    out = [cols[0], cols[2], log2ratio[cols[0]]]
                    print("\t".join(out), file=inactivated_out)

            elif cols[2] == "Escape":
                if cols[0] in log2ratio:
                    out = [cols[0], cols[2], log2ratio[cols[0]]]
                    print("\t".join(out), file=escape_out)

            elif cols[2] == "Variable_escape":
                if cols[0] in log2ratio:
                    out = [cols[0], cols[2], log2ratio[cols[0]]]
                    print("\t".join(out), file=variable_out)

