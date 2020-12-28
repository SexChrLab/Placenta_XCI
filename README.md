# Placenta_XCI
Understanding patterns of X chromosome inactivation in full term human placenta

### Processing placenta data
- Directory `placenta`
#### 01_process_dna
- Generate config files
    - `python generate_json_config_dna_females.py`: take in the input file `female_sample_ids.csv` and output the config file `process_dna_females_config.json` (for female placentas)
    - `python generate_json_config_dna_males.py`: take in the input file `male_sample_ids.csv` and output the config file `process_dna_males_config.json` (for male placenta)
    
- Snakemake files:
    - for mapping and genotype variants
    - `process_dna_females.snakefile` and `process_dna_males.snakefile`:
 
#### 02_run_asereadcounter
- Config file: `asereadcounter_config.json`
- Snakemake file: `asereadcounter.snakefile`
- Output directory: `02_run_asereadcounter/asereadcounter`

#### 03_analyze_ase
- Subdirectories: `scripts` and `results`
- Calculate (unphased) allele balance:
    - Python script `calc_allele_balance.py`
    - Config file: `analyze_ase_config.json`
    - Snakemake file: `analyze_ase.snakefile`
- Calculate median allele balance per individual:
    - Python script `calc_median_allele_balance_placenta_decidua.py` (see Bash script `run_calc_median_allele_balance_placenta_decidua.sh`)
    
#### 04_phasing
- Phasing strategy:
    - For each pair of placenta (site A and site B):
        + Subset to contain shared expressed variants
        + Using the site with more variants where allele balance is greater than 0.8
        + Generate a haplotype by adding all the biased allele together. If the allele balance is equal to 0.5, pick at random
        + Calculate allele balance using the phased data
- Steps:
  1. For each pair of placenta (site A and B), find shared variants between site A and site B:
    1. `python subset_paired_placentas_for_shared_variants.py chrX > chrX_summary_stats.txt`
    2. `python subset_paired_placentas_for_shared_variants.py chr8 > chr8_summary_stats.txt`
    2. Results are in directory `paired_placentas_shared_variants/`
  2. Run snakefile: `snakemake --snakefile phase.snakefile` to compute allele balance for phased data
  3. Concat for plotting Figure 2:
    ```
    cd 04_phasing/phased_allele_balance/
    cat *chrX*allele_balance_summary.tsv | grep -v sample_id | sort -n -r -k 3,3 > all_placenta_chrX_phased_allele_balance.tsv
    ```

#### 05_pca
- Contain files for generating the PCA

### Processing GTEx data
- Directory `gtex`
- In this directory, we are analyzing the ASEReadCounter counts from GTEx version 8.
## 01_download_data
1. Download the file `participant.tsv` from anvil project website. This file has information about the sample id
2. Download the file `sample.tsv` from anvil project website. This file has information about the rna id and which tissue
3. Obtain a list of individuals
  1. There are 979 individuals
  2. Run the python script `obtain_individuals_list.py`
4. Download using the file `download_asereadcounter_count.sh`
5. After downloading, I noticed that there are some files with this message inside: No such object: fc-secure-ff8156a3-ddf3-42e4-9211-0fd89da62108/GTEx_Analysis_2017-06-05_v8_ASE_WASP_chrX_raw_counts_by_subject/GTEX-1J8EW.v8.readcounts.chrX.txt.gz. We want to remove these individuals from further analyses. Therefore, we need to know which are these individuals.
  1. Run the python script: `python check_corrupted_files.py`
  2. The outfile is `failed_files.txt`. There are 147 individuals without ASEReadCounter results.
6. Generate a config file from the file `sample.tsv`: `python generate_config.py`
  1. Remove the individuals in the `failed_files.txt`
  2. Only keep the females
    1. Download: `wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt`

## 02_analyze_gtex_counts
- Snakemake file: `analyze_gtex_counts.snakefile`
1. Subset each downloaded count file for each tissue
  - Because each count file includes all of the tissues for an individual, I need to subset for each tissue
  - Use the python script `subset_gtex_counts.py`. See snakemake rule `subset_gtex_counts` (line 9)
2. calculate allele balance
  1.See snakemake rule `calc_allele_balance` (line 24)
3. Calculate median allele balance per tissue:
  1. See snakemake rule `calc_median_allele_balance_per_tissue` (line 36)
4. Find tissues where there are at least 10 samples per tissue: `python find_tissues_more_than_10_samples_per_tissue.py`

### Processing GTEx heart data
- Directory `heart`
## 1. Find subject ids with rnaseq data for both heart left ventricle and heart atrial appendage
- `python find_samples_with_2_hearts.py`

## 2. Calculate the proportion of skewed variants per sample
- `python calc_prop_variants_skewed_per_sample.py`

## 3. Employ a phasing strategy for heart in the same way as placenta
1. `python /scratch/tphung3/Placenta_XCI/heart/subset_paired_hearts_for_shared_variants.py chrX`. Results are in directory `paired_hearts_shared_variants/`
2. Use the snakemake file `phase.snakefile`. Results are in directory `phased_allele_balance/`
- `cat *chrX* | grep -v sample_id | awk '{print$1"\t"$3"\t"$2}' | sort -n -r -k 3,3 > all_heart_chrX_phased_allele_balance.tsv`

### Analysis at the gene level
- In this directory, I am analyzing genes that escape XCI for placenta and gtex tissues using the individuals that show skewed allele balance (median allele balance is greater than 0.8)
- Directory `gene_level`
    - Sub-directories: `gtex_counts` and `specific_gene_analysis`

#### wes_genotyping

#### gtex_counts
1. Find samples that are skewed in GTEX tissues, placenta, decidua females, and decidua males
  ```
  python /scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/scripts/find_samples_skewed.py
  ```
  - The config file is `/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/escape_genes_config.json`
  - There are 525 samples that are skewed in the GTEx
2. Tabulate how many individuals there are for each tissue that are highly skewed
  ```
  python /scratch/tphung3/PlacentaSexDiff/E_escape_genes/gtex_counts/scripts/tabulate_individuals.py
  ```
  - The result file is `tabulate_individuals.csv`
3. Convert skewed samples to bed file format
  1. Rule `convert_asereadcounter_to_bed_gtex` in `escape_genes.snakefile`
4. Use bedtools to find where on the genes the variants are
  1. Snakemake rule `bedtools_intersect`
  2. Remove duplicated: snakemake rule `find_unique_lines_after_bedtools`
  
##### Assign escape status by combining all of the gtex samples
1. Find genes that have at least one heterozygous and expressed variant across all skewed samples
  1. `python /scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/scripts/find_expressed_genes.py` produces the output file `/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/expressed_genes_all_samples.txt` that lists all of the genes with at least one heterozygous and expressed variant across all samples. There are 689 genes.
2. For each skewed individual in gtex, compute mean allele balance for each gene. For the placenta and decidua sample, I have already done this step here `/scratch/tphung3/Placenta_XCI/gene_level/wes_genotyping/asereadcounter_geneinfo/chrX`.
  - Use the Python script `compute_allele_balance_per_gene.py`
  - See snakemake rule `compute_allele_balance_per_gene_gtex`
3. Run `python make_allele_count_per_gene.py`
4. Add genes to config: `python add_genes_to_config.py`
5. See snakemake rule `plot_per_gene_allele_balance_compare_gtex_placenta_decidua`
6. Categorize genes into inactivated, escape, or variable escape for gtex, placenta, decidua females, and decidua males.
  1. Use Python script `categorize_genes.py`
    - This script categorizes the genes, remove NA, and also sort for plotting heatmaps

### Generate plots for manuscript
- Directory `manuscript_plots`
- Figure 2: `scripts/figure_2.R`
- Figure 3B: `scripts/figure_3B.R`