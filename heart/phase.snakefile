configfile: "/scratch/tphung3/Placenta_XCI/heart/subjects_with_two_hearts_config.json"

chromosomes = ["chrX"]

rule all:
    input:
        expand("phased_allele_balance/{subject}_{chr}_phased_allele_balance.tsv", subject=config["subjectids"], chr=chromosomes)

rule phase:
    input:
        "paired_hearts_shared_variants/{subject}_paired_hearts_shared_variants_{chr}.csv"
    output:
        "phased_allele_balance/{subject}_{chr}_phased_allele_balance.tsv"
    params:
        script = "/scratch/tphung3/Placenta_XCI/placenta/04_phasing/phase.py",
        subject = "{subject}",
        chromosome = "{chr}"
    shell:
        """
        python {params.script} --sample_id {params.subject} --chromosome {params.chromosome} --shared_variants {input} --outfile {output}
        """
