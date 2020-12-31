import os

configfile: "/scratch/tphung3/Placenta_XCI/placenta/01_process_dna/process_dna_females_config.json"

# conda environment: placenta

# Tool paths:
bwa_path = "bwa"
samtools_path = "samtools"

rule all:
    input:
        expand("bams/{sample}.GRCh38.default.sorted.bam.bai", sample=config["all_dna_samples"]),
        expand("/scratch/tphung3/Placenta_XCI/placenta/00_sexcheck/sexcheck_out/{sample}/{sample}_summary.tsv", sample=config["all_dna_samples"])

rule map:
    input:
        fq1 = lambda wildcards: os.path.join(config[wildcards.sample]["fq_path"], config[wildcards.sample]["fq_1"]),
        fq2 = lambda wildcards: os.path.join(config[wildcards.sample]["fq_path"], config[wildcards.sample]["fq_2"]),
        ref = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa",
        fai = "/data/CEM/shared/public_data/references/GENCODE/GRCh38.p12.genome/GRCh38.p12.genome.fa.fai"
    output:
        "bams/{sample}.GRCh38.default.sorted.bam"
    params:
        id = lambda wildcards: config[wildcards.sample]["ID"],
        sm = lambda wildcards: config[wildcards.sample]["SM"],
        lb = lambda wildcards: config[wildcards.sample]["LB"],
        pu = lambda wildcards: config[wildcards.sample]["PU"],
        pl = lambda wildcards: config[wildcards.sample]["PL"],
        bwa = bwa_path,
        samtools = samtools_path,
        threads = 4
    threads: 4
    priority: 100
    shell:
        " {params.bwa} mem -R "
        "'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
        "{input.ref} {input.fq1} {input.fq2}"
        "| {params.samtools} fixmate -O bam - - | {params.samtools} sort "
        "-O bam -o {output}"

rule index_bam:
    input:
        "bams/{sample}.GRCh38.default.sorted.bam"
    output:
        "bams/{sample}.GRCh38.default.sorted.bam.bai"
    params:
        samtools = samtools_path
    shell:
        "{params.samtools} index {input}"

rule sex_check:
    input:
        "bams/{sample}.GRCh38.default.sorted.bam"
    output:
        "/scratch/tphung3/Placenta_XCI/placenta/00_sexcheck/sexcheck_out/{sample}/{sample}_summary.tsv"
    params:
        sample_id = "{sample}",
        out_dir = "/scratch/tphung3/Placenta_XCI/placenta/00_sexcheck/sexcheck_out"
    shell:
        """
        python /home/tphung3/softwares/tanya_repos/SexInference/DNAseq/infer_sex_from_readsmapped.py --sample {params.sample_id} --bam_path {input} --out_dir {params.out_dir}
        """
