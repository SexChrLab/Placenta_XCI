#!/bin/bash
#SBATCH --job-name=escape_genes  # Job name
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=tphung3@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00

snakemake --snakefile escape_genes.snakefile -j 69 --cluster "sbatch --mem=12000 -c 4 -t 20:00:00" --rerun-incomplete
