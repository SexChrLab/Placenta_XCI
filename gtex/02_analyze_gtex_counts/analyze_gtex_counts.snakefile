# README
# Replace the paths on lines 30, 44, 71, etc...
import os

configfile: "/scratch/tphung3/Placenta_XCI/gtex/01_download_data/gtex_counts_analysis_config.json"

rule all:
    # autosomes
    input:
        expand("/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/autosomes_WASP_raw_counts_allele_balance/chr8/{tissue}_chr8_median_allele_balance.txt", tissue = config["tissues_ids"])
    input:
        expand("autosomes_WASP_raw_counts_allele_balance/chr8/{rna_id}_allele_balance.tsv", rna_id=config["all_rna_ids"])
    input:
        expand("autosomes_WASP_raw_counts/{rna_id}.v8.readcounts.autosomes.txt", rna_id=config["all_rna_ids"])
    # chrX
    input: #calculate median allele balance per tissue
        expand("/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/{tissue}_chrX_nonpars_median_allele_balance.txt", tissue = config["tissues_ids"])
    input: #calculate allele balance
        expand("chrX_WASP_raw_counts_allele_balance/{rna_id}_allele_balance.tsv", rna_id=config["all_rna_ids"])
    input: #subset for each rna_id
        expand("chrX_WASP_raw_counts/{rna_id}.v8.readcounts.chrX.txt", rna_id=config["all_rna_ids"])

# chrX
rule subset_gtex_counts_chrX:
    input:
    output:
        "chrX_WASP_raw_counts/{rna_id}.v8.readcounts.chrX.txt"
    params:
        rna_id = "{rna_id}",
        counts_dir = "/data/mwilsons/public_data/controlled_access/gtex/version8/chrX_WASP_raw_counts/",
        ending_str = ".v8.readcounts.chrX.txt.gz",
        chr = "chrX"
    shell:
        """
        python /scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/subset_gtex_counts.py --rna_id {params.rna_id} --counts_dir {params.counts_dir} --ending_str {params.ending_str} --chr {params.chr} --outfile {output}
        """

rule calc_allele_balance_chrX:
    input:
        "chrX_WASP_raw_counts/{rna_id}.v8.readcounts.chrX.txt"
    output:
        "chrX_WASP_raw_counts_allele_balance/{rna_id}_allele_balance.tsv"
    params:
        script = "/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/calc_allele_balance.py"
    shell:
        """
        python {params.script} --input {input} --output {output}
        """

rule calc_median_allele_balance_per_tissue_chrX:
    output:
        "/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/{tissue}_chrX_nonpars_median_allele_balance.txt"
    params:
        script = "/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/calc_median_allele_balance_per_tissue.py",
        tissue = "{tissue}",
        chr = "chrX",
        dir = "/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/",
        threshold = 10
    shell:
        """
        python {params.script} --tissue {params.tissue} --chromosome {params.chr} --files_dir {params.dir} --threshold {params.threshold} --include_pars no --out {output}
        """

# Autosomes
rule subset_gtex_counts_autosomes:
    input:
    output:
        "autosomes_WASP_raw_counts/{rna_id}.v8.readcounts.autosomes.txt"
    params:
        rna_id = "{rna_id}",
        counts_dir = "/data/mwilsons/public_data/controlled_access/gtex/version8/autosomes_WASP_raw_counts/",
        ending_str = ".v8.wasp_corrected.ase_table.tsv.gz",
        chr = "autosomes"
    shell:
        """
        python /scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/subset_gtex_counts.py --rna_id {params.rna_id} --counts_dir {params.counts_dir} --ending_str {params.ending_str} --chr {params.chr} --outfile {output}
        """

rule calc_allele_balance_autosomes:
    input:
        "autosomes_WASP_raw_counts/{rna_id}.v8.readcounts.autosomes.txt"
    output:
        "autosomes_WASP_raw_counts_allele_balance/chr8/{rna_id}_allele_balance.tsv"
    params:
        script = "/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/calc_allele_balance_gtex_counts_autosomes.py",
        chr = "chr8"
    shell:
        """
        python {params.script} --input {input} --chr {params.chr} --output {output}
        """

rule calc_median_allele_balance_per_tissue_autosomes:
    output:
        "/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/autosomes_WASP_raw_counts_allele_balance/chr8/{tissue}_chr8_median_allele_balance.txt"
    params:
        script = "/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/calc_median_allele_balance_per_tissue.py",
        tissue = "{tissue}",
        chr = "chr8",
        dir = "/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/autosomes_WASP_raw_counts_allele_balance/chr8",
        threshold = 10
    shell:
        """
        python {params.script} --tissue {params.tissue} --chromosome {params.chr} --files_dir {params.dir} --threshold {params.threshold} --include_pars yes --out {output}
        """
