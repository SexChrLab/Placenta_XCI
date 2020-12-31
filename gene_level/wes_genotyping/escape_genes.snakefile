import os

configfile: "/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/escape_genes_config.json"

rule all:
    input:
        expand("/scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/asereadcounter_geneinfo/chrX/{sample}_chrX_allele_balance_geneinfo_unique.tsv", sample=config["full_term_placenta"])

rule convert_asereadcounter_to_bed_placenta:
    input:
        "/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chrX/{sample}_chrX_allele_balance.tsv"
    output:
        "/scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/bed/chrX/{sample}_chrX_allele_balance.bed"
    shell:
        """
        python /scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/convert_asereadcounter_to_bed.py --input_asereadcounter {input} --output_bed {output}
        """

rule bedtools_intersect:
    input:
        "/scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/bed/chrX/{sample}_chrX_allele_balance.bed"
    output:
        "/scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/asereadcounter_geneinfo/chrX/{sample}_chrX_allele_balance_geneinfo.tsv"
    shell:
        """
        bedtools intersect -a /scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/gtf_bed/gencode.v29.annotation.chrX.bed -b {input} -wa -wb > {output}
        """

rule find_unique_lines_after_bedtools:
    input:
        "/scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/asereadcounter_geneinfo/chrX/{sample}_chrX_allele_balance_geneinfo.tsv"
    output:
        "/scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/asereadcounter_geneinfo/chrX/{sample}_chrX_allele_balance_geneinfo_unique.tsv"
    shell:
        """
        python /scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/find_unique_lines_after_bedtools.py --infile {input} --outfile {output} --index 5
        """
