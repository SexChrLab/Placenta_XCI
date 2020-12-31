# In this script, we want to calculate: for each rna, the proportion of variants that are skewed
import json
import pandas as pd

threshold = 10

# Placenta
placenta_out = open("/scratch/tphung3/Placenta_XCI/heart/placenta_skewed_prop.txt", "w")
with open("/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/scripts/analyze_ase_config.json") as json_file:
    data = json.load(json_file)
    for sampleid in data["dna_rna_placenta"]:
        # site 1
        site_1_id = sampleid + "_" + data["dna_rna_placenta"][sampleid][0]
        site_1 = pd.read_csv("/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chrX/" + site_1_id + "_chrX_allele_balance.tsv", sep="\t")
        site_1_subset = site_1[site_1['total_count'] > threshold]
        site_1_subset_nrow = len(site_1_subset.index)

        site_1_subset_skewed = site_1_subset[site_1_subset['allele_balance'] > 0.8]
        site_1_subset_skewed_nrow = len(site_1_subset_skewed.index)

        site_1_skewed_prop = site_1_subset_skewed_nrow/site_1_subset_nrow

        # site 2
        site_2_id = sampleid + "_" + data["dna_rna_placenta"][sampleid][1]
        site_2 = pd.read_csv("/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chrX/" + site_2_id + "_chrX_allele_balance.tsv", sep="\t")
        site_2_subset = site_2[site_2['total_count'] > threshold]
        site_2_subset_nrow = len(site_2_subset.index)

        site_2_subset_skewed = site_2_subset[site_2_subset['allele_balance'] > 0.8]
        site_2_subset_skewed_nrow = len(site_2_subset_skewed.index)

        site_2_skewed_prop = site_2_subset_skewed_nrow/site_2_subset_nrow

        if site_1_skewed_prop > site_2_skewed_prop:
            out = [sampleid, str(site_1_skewed_prop)]
            print("\t".join(out), file=placenta_out)
        else:
            out = [sampleid, str(site_2_skewed_prop)]
            print("\t".join(out), file=placenta_out)

# Heart
heart_out = open("/scratch/tphung3/Placenta_XCI/heart/heart_skewed_prop.txt", "w")
with open("/scratch/tphung3/Placenta_XCI/heart/subjects_with_two_hearts.txt", "r") as f:
    for line in f:
        if not line.startswith("subjectid"):
            items = line.rstrip("\n").split("\t")

            # heart_left_ventricle
            heart_left_ventricle = pd.read_csv("/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/" + items[1] + "_allele_balance.tsv", sep="\t")
            heart_left_ventricle_subset = heart_left_ventricle[heart_left_ventricle['total_count'] > threshold]
            heart_left_ventricle_subset_nrow = len(heart_left_ventricle_subset.index)

            heart_left_ventricle_subset_skewed = heart_left_ventricle_subset[heart_left_ventricle_subset['allele_balance'] > 0.8]
            heart_left_ventricle_subset_skewed_nrow = len(heart_left_ventricle_subset_skewed.index)

            heart_left_ventricle_skewed_prop = heart_left_ventricle_subset_skewed_nrow/heart_left_ventricle_subset_nrow

            # heart_atrial_appendage
            heart_atrial_appendage = pd.read_csv("/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/" + items[2] + "_allele_balance.tsv", sep="\t")
            heart_atrial_appendage_subset = heart_atrial_appendage[heart_atrial_appendage['total_count'] > threshold]
            heart_atrial_appendage_subset_nrow = len(heart_atrial_appendage_subset.index)

            heart_atrial_appendage_subset_skewed = heart_atrial_appendage_subset[heart_atrial_appendage_subset['allele_balance'] > 0.8]
            heart_atrial_appendage_subset_skewed_nrow = len(heart_atrial_appendage_subset_skewed.index)

            heart_atrial_appendage_skewed_prop = heart_atrial_appendage_subset_skewed_nrow/heart_atrial_appendage_subset_nrow

            if heart_left_ventricle_skewed_prop > heart_atrial_appendage_skewed_prop:
                out = [items[1], str(heart_left_ventricle_skewed_prop)]
                print("\t".join(out), file=heart_out)
            else:
                out = [items[1], str(heart_atrial_appendage_skewed_prop)]
                print("\t".join(out), file=heart_out)




