# In this script, we are finding individuals with skewed allele balance
import json
import os

data = {}
data['all_samples'] = []
data['gtex_samples'] = []
data['full_term_placenta'] = []
data['decidua_females'] = []
data['decidua_males'] = []
data['tissues_ids'] = []
data['tissues_ids_for_gene_analysis'] = []

# GTEX tissues
gtex_tissues = ["Thyroid", "Uterus", "Pancreas", "Artery_Aorta", "Spleen", "Small_Intestine_Terminal_Ileum", "Skin_Not_Sun_Exposed_Suprapubic", "Pituitary", "Esophagus_Mucosa", "Brain_Cortex", "Artery_Tibial", "Artery_Coronary", "Kidney_Cortex", "Nerve_Tibial", "Brain_Cerebellum", "Lung", "Adrenal_Gland", "Heart_Left_Ventricle", "Brain_Caudate_basal_ganglia", "Esophagus_Muscularis", "Esophagus_Gastroesophageal_Junction", "Brain_Cerebellar_Hemisphere", "Brain_Spinal_cord_cervical_c-1", "Brain_Frontal_Cortex_BA9", "Brain_Hypothalamus", "Minor_Salivary_Gland", "Whole_Blood", "Brain_Substantia_nigra", "Liver", "Colon_Transverse", "Heart_Atrial_Appendage", "Brain_Amygdala", "Stomach", "Ovary", "Adipose_Visceral_Omentum", "Skin_Sun_Exposed_Lower_leg", "Colon_Sigmoid", "Brain_Anterior_cingulate_cortex_BA24", "Adipose_Subcutaneous", "Brain_Hippocampus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Vagina", "Muscle_Skeletal", "Breast_Mammary_Tissue"]

for tissue in gtex_tissues:
    skewed_samples = []
    with open('/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/' + tissue + '_chrX_nonpars_median_allele_balance.txt', 'r') as f:
        for line in f:
            if not line.startswith('tissue'):
                items = line.rstrip('\n').split('\t')
                if float(items[6]) > 0.8:
                    id = items[2]
                    skewed_samples.append(id)
                    data['all_samples'].append(id)
                    data['gtex_samples'].append(id)
    if len(skewed_samples) != 0:
        data[tissue] = skewed_samples
        data['tissues_ids'].append(tissue)
    if len(skewed_samples) >= 10:
        data['tissues_ids_for_gene_analysis'].append(tissue)

# print (data['tissues_ids_for_gene_analysis'])
# print (len(data['tissues_ids_for_gene_analysis']))
print(len(data['gtex_samples']))

placenta_skewed_samples = []
with open('/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chrX/placenta_chrX_nonpars_median_allele_balance.txt', 'r') as f:
    for line in f:
        if not line.startswith('tissue'):
            items = line.rstrip('\n').split('\t')
            if items[2] != "OBG0175": #remove )BG0175 samples
                if float(items[7]) > 0.8:
                    id = items[2] + '_' + items[3]
                    placenta_skewed_samples.append(id)
                    data['all_samples'].append(id)
                    data['full_term_placenta'].append(id)

# Comment out the deciduas
# decidua_females_skewed_samples = []
# with open('/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chrX/decidua_females_chrX_nonpars_median_allele_balance.txt', 'r') as f:
#     for line in f:
#         if not line.startswith('tissue'):
#             items = line.rstrip('\n').split('\t')
#             if float(items[7]) > 0.8:
#                 id = items[2] + '_' + items[3]
#                 decidua_females_skewed_samples.append(id)
#                 data['all_samples'].append(id)
#                 data['decidua_females'].append(id)
#
# decidua_males_skewed_samples = []
# with open('/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chrX/decidua_males_chrX_nonpars_median_allele_balance.txt', 'r') as f:
#     for line in f:
#         if not line.startswith('tissue'):
#             items = line.rstrip('\n').split('\t')
#             if float(items[7]) > 0.8:
#                 id = items[2] + '_' + items[3]
#                 decidua_males_skewed_samples.append(id)
#                 data['all_samples'].append(id)
#                 data['decidua_males'].append(id)

with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/escape_genes_config.json', 'w') as outfile:
    json.dump(data, outfile)

# Count how many placenta samples are skewed
# with open('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/escape_genes_config.json') as json_file:
#     data = json.load(json_file)
#     print(len(data['full_term_placenta']))

