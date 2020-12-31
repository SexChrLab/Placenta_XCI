configfile: "/data/CEM/wilsonlab/projects/placenta/Placenta_XCI/placenta/03_analyze_ase/scripts/analyze_ase_config.json"

chromosomes = ["chr8", "chrX"]

rule all:
    input:
        expand("phased_allele_balance/{sample}_{chr}_phased_allele_balance_summary.tsv", sample=config["dna_rna_placenta"], chr=chromosomes),
        expand("phased_allele_balance/{sample}_{chr}_phased_allele_balance_data.tsv", sample=config["dna_rna_placenta"], chr=chromosomes)

rule phase:
    input:
        "/data/CEM/wilsonlab/projects/placenta/Placenta_XCI/placenta/04_phasing/paired_placentas_shared_variants/{sample}_paired_placentas_shared_variants_{chr}.csv"
    output:
        summary = "phased_allele_balance/{sample}_{chr}_phased_allele_balance_summary.tsv",
        phased_data = "phased_allele_balance/{sample}_{chr}_phased_allele_balance_data.tsv"
    params:
        script = "/data/CEM/wilsonlab/projects/placenta/Placenta_XCI/placenta/04_phasing/phase.py",
        sample = "{sample}",
        chromosome = "{chr}"
    shell:
        """
        python {params.script} --sample_id {params.sample} --chromosome {params.chromosome} --shared_variants {input} --out_summary {output.summary} --out_phased_data {output.phased_data}
        """
