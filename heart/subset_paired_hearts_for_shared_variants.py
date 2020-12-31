# This is the first script in phasing where we are subsetting the paired placentas for shared variants only
import pandas as pd
import json
import sys

chr = sys.argv[1]

totalCount_threshold = 10 #Change this if needs be

with open("/scratch/tphung3/Placenta_XCI/heart/subjects_with_two_hearts.txt") as f:
    for line in f:
        if not line.startswith("subjectid"):
            items = line.rstrip("\n").split("\t")
            heart_left_ventricle = pd.read_csv(
                "/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/" +
                items[1] + "_allele_balance.tsv", sep="\t")
            heart_atrial_appendage = pd.read_csv(
                "/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/" +
                items[2] + "_allele_balance.tsv", sep="\t")

            # Find shared variants
            shared_variants = pd.merge(heart_left_ventricle, heart_atrial_appendage, on='position', how='inner')

            # Filter based on total_count threshold
            shared_variants_threshold = shared_variants[(shared_variants['total_count_x']>totalCount_threshold) & (shared_variants['total_count_y']>totalCount_threshold)]

            shared_variants_threshold.to_csv("/scratch/tphung3/Placenta_XCI/heart/paired_hearts_shared_variants/" + items[0] + "_paired_hearts_shared_variants_" + chr + ".csv", index=False)