import os

configfile: "/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/escape_genes_config.json"

rule all:
    input: #plotting
        expand("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/plots/{gene}_allele_balance_compare_gtex_placenta_decidua_violin.png", gene=config["genes"])
    input: #compute allele balance per gene for gtex samples
        expand("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/asereadcounter_geneinfo/chrX/{sample}_chrX_allele_balance_per_gene.tsv", sample=config["gtex_samples"])
    input: #bedtools and find unique
        expand("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/asereadcounter_geneinfo/chrX/{sample}_chrX_allele_balance_geneinfo_unique.tsv", sample=config["gtex_samples"])
    input: #convert gtex asereadcounter to bed file format
        expand("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/bed/chrX/{sample}_chrX_allele_balance.bed", sample=config["gtex_samples"])


rule convert_asereadcounter_to_bed_gtex:
    input:
        "/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/{sample}_allele_balance.tsv"
    output:
        "/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/bed/chrX/{sample}_chrX_allele_balance.bed"
    shell:
        """
        python /scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/scripts/convert_asereadcounter_to_bed.py --input_asereadcounter {input} --output_bed {output}
        """

rule bedtools_intersect:
    input:
        "/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/bed/chrX/{sample}_chrX_allele_balance.bed"
    output:
        "/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/asereadcounter_geneinfo/chrX/{sample}_chrX_allele_balance_geneinfo.tsv"
    shell:
        """
        bedtools intersect -a /scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/misc/gtf_bed/gencode.v29.annotation.chrX.bed -b {input} -wa -wb > {output}
        """

rule find_unique_lines_after_bedtools:
    input:
        "/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/asereadcounter_geneinfo/chrX/{sample}_chrX_allele_balance_geneinfo.tsv"
    output:
        "/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/asereadcounter_geneinfo/chrX/{sample}_chrX_allele_balance_geneinfo_unique.tsv"
    shell:
        """
        python /scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/scripts/find_unique_lines_after_bedtools.py --infile {input} --outfile {output} --index 5
        """

rule compute_allele_balance_per_gene_gtex:
    input:
        "/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/asereadcounter_geneinfo/chrX/{sample}_chrX_allele_balance_geneinfo_unique.tsv"
    output:
        "/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/asereadcounter_geneinfo/chrX/{sample}_chrX_allele_balance_per_gene.tsv"
    shell:
        """
        python /scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/scripts/compute_allele_balance_per_gene.py --input_asereadcounter {input} --output {output} --threshold 10
        """

rule plot_per_gene_allele_balance_compare_gtex_placenta_decidua:
    input:
    output:
        "/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/plots/{gene}_allele_balance_compare_gtex_placenta_decidua_violin.png"
    params:
        gene = "{gene}"
    shell:
        """
        Rscript /scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/plot_per_gene_allele_balance_compare_gtex_placenta_violin.R {params.gene}
        """
