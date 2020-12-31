# In this script, we want to output:
# for each individual, what are the tissues thare are skewed
from collections import defaultdict
gtex_tissues = ['Adipose_Subcutaneous', 'Adipose_Visceral_Omentum', 'Adrenal_Gland', 'Artery_Aorta', 'Artery_Coronary', 'Artery_Tibial', 'Brain_Amygdala', 'Brain_Anterior_cingulate_cortex_BA24', 'Brain_Caudate_basal_ganglia', 'Brain_Cerebellar_Hemisphere', 'Brain_Cerebellum', 'Brain_Cortex', 'Brain_Frontal_Cortex_BA9', 'Brain_Hippocampus', 'Brain_Hypothalamus', 'Brain_Nucleus_accumbens_basal_ganglia', 'Brain_Putamen_basal_ganglia', 'Brain_Spinal_cord_cervical_c-1', 'Brain_Substantia_nigra', 'Breast_Mammary_Tissue', 'Colon_Sigmoid', 'Colon_Transverse', 'Esophagus_Gastroesophageal_Junction', 'Esophagus_Mucosa', 'Esophagus_Muscularis', 'Heart_Atrial_Appendage', 'Heart_Left_Ventricle', 'Kidney_Cortex', 'Liver', 'Lung', 'Minor_Salivary_Gland', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Skin_Not_Sun_Exposed_Suprapubic', 'Skin_Sun_Exposed_Lower_leg', 'Small_Intestine_Terminal_Ileum', 'Spleen', 'Stomach', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood']

skewed_rnaid_tissue = {}

for tissue in gtex_tissues:
    with open('/scratch/tphung3/PlacentaSexDiff/F_gtex_counts/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/' + tissue + '_chrX_nonpars_median_allele_balance.txt', 'r') as f:
        for line in f:
            if not line.startswith('tissue'):
                items = line.rstrip('\n').split('\t')
                if float(items[6]) > 0.8:
                    skewed_rnaid_tissue[items[2]] = tissue

skewed_sampleid_tissue = defaultdict(list)
for skewed_rnaid in skewed_rnaid_tissue:
    items = skewed_rnaid.split('-')
    sampleid = items[0] + '-' + items[1]
    skewed_sampleid_tissue[sampleid].append(skewed_rnaid_tissue[skewed_rnaid])

for k, v in skewed_sampleid_tissue.items():
    if len(v) > 2:
        out = [k]
        for i in sorted(v):
            out.append(i)
        print(','.join(out))
