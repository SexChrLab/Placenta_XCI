# chr8
setwd('/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/autosomes_WASP_raw_counts_allele_balance/chr8/')
chr8_Adipose_Subcutaneous = read.csv('Adipose_Subcutaneous_chr8_median_allele_balance.txt', sep = '\t')
chr8_Adipose_Visceral_Omentum = read.csv('Adipose_Visceral_Omentum_chr8_median_allele_balance.txt', sep = '\t')
chr8_Adrenal_Gland = read.csv('Adrenal_Gland_chr8_median_allele_balance.txt', sep = '\t')
chr8_Artery_Aorta = read.csv('Artery_Aorta_chr8_median_allele_balance.txt', sep = '\t')
chr8_Artery_Coronary = read.csv('Artery_Coronary_chr8_median_allele_balance.txt', sep = '\t')
chr8_Artery_Tibial = read.csv('Artery_Tibial_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Amygdala = read.csv('Brain_Amygdala_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Anterior_cingulate_cortex_BA24 = read.csv('Brain_Anterior_cingulate_cortex_BA24_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Caudate_basal_ganglia = read.csv('Brain_Caudate_basal_ganglia_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Cerebellar_Hemisphere = read.csv('Brain_Cerebellar_Hemisphere_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Cerebellum = read.csv('Brain_Cerebellum_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Cortex = read.csv('Brain_Cortex_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Frontal_Cortex_BA9 = read.csv('Brain_Frontal_Cortex_BA9_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Hippocampus = read.csv('Brain_Hippocampus_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Hypothalamus = read.csv('Brain_Hypothalamus_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Nucleus_accumbens_basal_ganglia = read.csv('Brain_Nucleus_accumbens_basal_ganglia_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Putamen_basal_ganglia = read.csv('Brain_Putamen_basal_ganglia_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Spinal_cord_cervical = read.csv('Brain_Spinal_cord_cervical_c-1_chr8_median_allele_balance.txt', sep = '\t')
chr8_Brain_Substantia_nigra = read.csv('Brain_Substantia_nigra_chr8_median_allele_balance.txt', sep = '\t')
chr8_Breast_Mammary_Tissue = read.csv('Breast_Mammary_Tissue_chr8_median_allele_balance.txt', sep = '\t')
chr8_Colon_Sigmoid = read.csv('Colon_Sigmoid_chr8_median_allele_balance.txt', sep = '\t')
chr8_Colon_Transverse = read.csv('Colon_Transverse_chr8_median_allele_balance.txt', sep = '\t')
chr8_Esophagus_Gastroesophageal_Junction = read.csv('Esophagus_Gastroesophageal_Junction_chr8_median_allele_balance.txt', sep = '\t')
chr8_Esophagus_Mucosa = read.csv('Esophagus_Mucosa_chr8_median_allele_balance.txt', sep = '\t')
chr8_Esophagus_Muscularis = read.csv('Esophagus_Muscularis_chr8_median_allele_balance.txt', sep = '\t')
chr8_Heart_Atrial_Appendage = read.csv('Heart_Atrial_Appendage_chr8_median_allele_balance.txt', sep = '\t')
chr8_Heart_Left_Ventricle = read.csv('Heart_Left_Ventricle_chr8_median_allele_balance.txt', sep = '\t')
chr8_Kidney_Cortex = read.csv('Kidney_Cortex_chr8_median_allele_balance.txt', sep = '\t')
chr8_Liver = read.csv('Liver_chr8_median_allele_balance.txt', sep = '\t')
chr8_Lung = read.csv('Lung_chr8_median_allele_balance.txt', sep = '\t')
chr8_Minor_Salivary_Gland = read.csv('Minor_Salivary_Gland_chr8_median_allele_balance.txt', sep = '\t')
chr8_Muscle_Skeletal = read.csv('Muscle_Skeletal_chr8_median_allele_balance.txt', sep = '\t')
chr8_Nerve_Tibial = read.csv('Nerve_Tibial_chr8_median_allele_balance.txt', sep = '\t')
chr8_Ovary = read.csv('Ovary_chr8_median_allele_balance.txt', sep = '\t')
chr8_Pancreas = read.csv('Pancreas_chr8_median_allele_balance.txt', sep = '\t')
chr8_Pituitary = read.csv('Pituitary_chr8_median_allele_balance.txt', sep = '\t')
chr8_Skin_Not_Sun_Exposed_Suprapubic = read.csv('Skin_Not_Sun_Exposed_Suprapubic_chr8_median_allele_balance.txt', sep = '\t')
chr8_Skin_Sun_Exposed_Lower_leg = read.csv('Skin_Sun_Exposed_Lower_leg_chr8_median_allele_balance.txt', sep = '\t')
chr8_Small_Intestine_Terminal_Ileum = read.csv('Small_Intestine_Terminal_Ileum_chr8_median_allele_balance.txt', sep = '\t')
chr8_Spleen = read.csv('Spleen_chr8_median_allele_balance.txt', sep = '\t')
chr8_Stomach = read.csv('Stomach_chr8_median_allele_balance.txt', sep = '\t')
chr8_Thyroid = read.csv('Thyroid_chr8_median_allele_balance.txt', sep = '\t')
chr8_Uterus = read.csv('Uterus_chr8_median_allele_balance.txt', sep = '\t')
chr8_Vagina = read.csv('Vagina_chr8_median_allele_balance.txt', sep = '\t')
chr8_Whole_Blood = read.csv('Whole_Blood_chr8_median_allele_balance.txt', sep = '\t')


# Placenta tissues
chr8_placenta = read.csv('/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chr8/placenta_chr8_median_allele_balance.txt', sep = '\t')
chr8_placenta_rm = subset(chr8_placenta, chr8_placenta$dna_id!="OBG0175")

Adipose_Subcutaneous_lab = rep("Adipose Subcutaneous", nrow(chr8_Adipose_Subcutaneous))
Adipose_Visceral_Omentum_lab = rep("Adipose Visceral Omentum", nrow(chr8_Adipose_Visceral_Omentum))
Adrenal_Gland_lab = rep("Adrenal Gland", nrow(chr8_Adrenal_Gland))
Artery_Aorta_lab = rep("Artery Aorta", nrow(chr8_Artery_Aorta))
Artery_Coronary_lab = rep("Artery Coronary", nrow(chr8_Artery_Coronary))
Artery_Tibial_lab = rep("Artery Tibial", nrow(chr8_Artery_Tibial))
Brain_Amygdala_lab = rep("Brain Amygdala", nrow(chr8_Brain_Amygdala))
Brain_Anterior_cingulate_cortex_BA24_lab = rep("Brain Anterior cingulate cortex BA24", nrow(chr8_Brain_Anterior_cingulate_cortex_BA24))
Brain_Caudate_basal_ganglia_lab = rep("Brain Caudate basal ganglia", nrow(chr8_Brain_Caudate_basal_ganglia))
Brain_Cerebellar_Hemisphere_lab = rep("Brain Cerebellar Hemisphere", nrow(chr8_Brain_Cerebellar_Hemisphere))
Brain_Cerebellum_lab = rep("Brain Cerebellum", nrow(chr8_Brain_Cerebellum))
Brain_Cortex_lab = rep("Brain Cortex", nrow(chr8_Brain_Cortex))
Brain_Frontal_Cortex_BA9_lab = rep("Brain Frontal Cortex BA9", nrow(chr8_Brain_Frontal_Cortex_BA9))
Brain_Hippocampus_lab = rep("Brain Hippocampus", nrow(chr8_Brain_Hippocampus))
Brain_Hypothalamus_lab = rep("Brain Hypothalamus", nrow(chr8_Brain_Hypothalamus))
Brain_Nucleus_accumbens_basal_ganglia_lab = rep("Brain Nucleus accumbens basal ganglia", nrow(chr8_Brain_Nucleus_accumbens_basal_ganglia))
Brain_Putamen_basal_ganglia_lab = rep("Brain Putamen basal ganglia", nrow(chr8_Brain_Putamen_basal_ganglia))
Brain_Spinal_cord_cervical_lab = rep("Brain Spinal cord cervical", nrow(chr8_Brain_Spinal_cord_cervical))
Brain_Substantia_nigra_lab = rep("Brain Substantia nigra", nrow(chr8_Brain_Substantia_nigra))
Breast_Mammary_Tissue_lab = rep("Breast Mammary Tissue", nrow(chr8_Breast_Mammary_Tissue))
Colon_Sigmoid_lab = rep("Colon Sigmoid", nrow(chr8_Colon_Sigmoid))
Colon_Transverse_lab = rep("Colon Transverse", nrow(chr8_Colon_Transverse))
Esophagus_Gastroesophageal_Junction_lab = rep("Esophagus Gastroesophageal Junction", nrow(chr8_Esophagus_Gastroesophageal_Junction))
Esophagus_Mucosa_lab = rep("Esophagus Mucosa", nrow(chr8_Esophagus_Mucosa))
Esophagus_Muscularis_lab = rep("Esophagus Muscularis", nrow(chr8_Esophagus_Muscularis))
Heart_Atrial_Appendage_lab = rep("Heart Atrial Appendage", nrow(chr8_Heart_Atrial_Appendage))
Heart_Left_Ventricle_lab = rep("Heart Left Ventricle", nrow(chr8_Heart_Left_Ventricle))
Kidney_Cortex_lab = rep("Kidney Cortex", nrow(chr8_Kidney_Cortex))
Liver_lab = rep("Liver", nrow(chr8_Liver))
Lung_lab = rep("Lung", nrow(chr8_Lung))
Minor_Salivary_Gland_lab = rep("Minor Salivary Gland", nrow(chr8_Minor_Salivary_Gland))
Muscle_Skeletal_lab = rep("Muscle Skeletal", nrow(chr8_Muscle_Skeletal))
Nerve_Tibial_lab = rep("Nerve Tibial", nrow(chr8_Nerve_Tibial))
Ovary_lab = rep("Ovary", nrow(chr8_Ovary))
Pancreas_lab = rep("Pancreas", nrow(chr8_Pancreas))
Pituitary_lab = rep("Pituitary", nrow(chr8_Pituitary))
Skin_Not_Sun_Exposed_Suprapubic_lab = rep("Skin Not Sun Exposed Suprapubic", nrow(chr8_Skin_Not_Sun_Exposed_Suprapubic))
Skin_Sun_Exposed_Lower_leg_lab = rep("Skin Sun Exposed Lower leg", nrow(chr8_Skin_Sun_Exposed_Lower_leg))
Small_Intestine_Terminal_Ileum_lab = rep("Small Intestine Terminal Ileum", nrow(chr8_Small_Intestine_Terminal_Ileum))
Spleen_lab = rep("Spleen", nrow(chr8_Spleen))
Stomach_lab = rep("Stomach", nrow(chr8_Stomach))
Thyroid_lab = rep("Thyroid", nrow(chr8_Thyroid))
Uterus_lab = rep("Uterus", nrow(chr8_Uterus))
Vagina_lab = rep("Vagina", nrow(chr8_Vagina))
Whole_Blood_lab = rep("Whole Blood", nrow(chr8_Whole_Blood))

placenta_lab = rep('Placenta', 58)

n_gtex = nrow(chr8_Thyroid)+ nrow(chr8_Uterus)+ nrow(chr8_Pancreas)+ nrow(chr8_Artery_Aorta)+ nrow(chr8_Spleen)+ nrow(chr8_Small_Intestine_Terminal_Ileum)+ nrow(chr8_Skin_Not_Sun_Exposed_Suprapubic)+ nrow(chr8_Pituitary)+ nrow(chr8_Esophagus_Mucosa)+ nrow(chr8_Brain_Cortex)+ nrow(chr8_Artery_Tibial)+ nrow(chr8_Artery_Coronary)+ nrow(chr8_Kidney_Cortex)+ nrow(chr8_Nerve_Tibial)+ nrow(chr8_Brain_Cerebellum)+ nrow(chr8_Lung)+ nrow(chr8_Adrenal_Gland)+ nrow(chr8_Heart_Left_Ventricle)+ nrow(chr8_Brain_Caudate_basal_ganglia)+ nrow(chr8_Esophagus_Muscularis)+ nrow(chr8_Esophagus_Gastroesophageal_Junction)+ nrow(chr8_Brain_Cerebellar_Hemisphere)+ nrow(chr8_Brain_Spinal_cord_cervical)+ nrow(chr8_Breast_Mammary_Tissue)+ nrow(chr8_Brain_Frontal_Cortex_BA9)+ nrow(chr8_Brain_Hypothalamus)+ nrow(chr8_Minor_Salivary_Gland)+ nrow(chr8_Whole_Blood)+ nrow(chr8_Brain_Substantia_nigra)+ nrow(chr8_Liver)+ nrow(chr8_Colon_Transverse)+ nrow(chr8_Heart_Atrial_Appendage)+ nrow(chr8_Brain_Amygdala)+ nrow(chr8_Stomach)+ nrow(chr8_Ovary)+ nrow(chr8_Adipose_Visceral_Omentum)+ nrow(chr8_Skin_Sun_Exposed_Lower_leg)+ nrow(chr8_Colon_Sigmoid)+ nrow(chr8_Brain_Anterior_cingulate_cortex_BA24)+ nrow(chr8_Adipose_Subcutaneous)+ nrow(chr8_Brain_Hippocampus)+ nrow(chr8_Brain_Nucleus_accumbens_basal_ganglia)+ nrow(chr8_Brain_Putamen_basal_ganglia)+ nrow(chr8_Vagina)+ nrow(chr8_Muscle_Skeletal)

gtex_lab = rep('GTEX', n_gtex)

chr8_data = data.frame(median = c(chr8_Adipose_Subcutaneous$median_allele_balance_subset,
                                  chr8_Adipose_Visceral_Omentum$median_allele_balance_subset,
                                  chr8_Adrenal_Gland$median_allele_balance_subset,
                                  chr8_Artery_Aorta$median_allele_balance_subset,
                                  chr8_Artery_Coronary$median_allele_balance_subset,
                                  chr8_Artery_Tibial$median_allele_balance_subset,
                                  chr8_Brain_Amygdala$median_allele_balance_subset,
                                  chr8_Brain_Anterior_cingulate_cortex_BA24$median_allele_balance_subset,
                                  chr8_Brain_Caudate_basal_ganglia$median_allele_balance_subset,
                                  chr8_Brain_Cerebellar_Hemisphere$median_allele_balance_subset,
                                  chr8_Brain_Cerebellum$median_allele_balance_subset,
                                  chr8_Brain_Cortex$median_allele_balance_subset,
                                  chr8_Brain_Frontal_Cortex_BA9$median_allele_balance_subset,
                                  chr8_Brain_Hippocampus$median_allele_balance_subset,
                                  chr8_Brain_Hypothalamus$median_allele_balance_subset,
                                  chr8_Brain_Nucleus_accumbens_basal_ganglia$median_allele_balance_subset,
                                  chr8_Brain_Putamen_basal_ganglia$median_allele_balance_subset,
                                  chr8_Brain_Spinal_cord_cervical$median_allele_balance_subset,
                                  chr8_Brain_Substantia_nigra$median_allele_balance_subset,
                                  chr8_Breast_Mammary_Tissue$median_allele_balance_subset,
                                  chr8_Colon_Sigmoid$median_allele_balance_subset,
                                  chr8_Colon_Transverse$median_allele_balance_subset,
                                  chr8_Esophagus_Gastroesophageal_Junction$median_allele_balance_subset,
                                  chr8_Esophagus_Mucosa$median_allele_balance_subset,
                                  chr8_Esophagus_Muscularis$median_allele_balance_subset,
                                  chr8_Heart_Atrial_Appendage$median_allele_balance_subset,
                                  chr8_Heart_Left_Ventricle$median_allele_balance_subset,
                                  chr8_Kidney_Cortex$median_allele_balance_subset,
                                  chr8_Liver$median_allele_balance_subset,
                                  chr8_Lung$median_allele_balance_subset,
                                  chr8_Minor_Salivary_Gland$median_allele_balance_subset,
                                  chr8_Muscle_Skeletal$median_allele_balance_subset,
                                  chr8_Nerve_Tibial$median_allele_balance_subset,
                                  chr8_Ovary$median_allele_balance_subset,
                                  chr8_Pancreas$median_allele_balance_subset,
                                  chr8_Pituitary$median_allele_balance_subset,
                                  chr8_Skin_Not_Sun_Exposed_Suprapubic$median_allele_balance_subset,
                                  chr8_Skin_Sun_Exposed_Lower_leg$median_allele_balance_subset,
                                  chr8_Small_Intestine_Terminal_Ileum$median_allele_balance_subset,
                                  chr8_Spleen$median_allele_balance_subset,
                                  chr8_Stomach$median_allele_balance_subset,
                                  chr8_Thyroid$median_allele_balance_subset,
                                  chr8_Uterus$median_allele_balance_subset,
                                  chr8_Vagina$median_allele_balance_subset,
                                  chr8_Whole_Blood$median_allele_balance_subset,
                                  chr8_placenta_rm$median_allele_balance_subset),
                       tissues = c(Adipose_Subcutaneous_lab,
                                   Adipose_Visceral_Omentum_lab,
                                   Adrenal_Gland_lab,
                                   Artery_Aorta_lab,
                                   Artery_Coronary_lab,
                                   Artery_Tibial_lab,
                                   Brain_Amygdala_lab,
                                   Brain_Anterior_cingulate_cortex_BA24_lab,
                                   Brain_Caudate_basal_ganglia_lab,
                                   Brain_Cerebellar_Hemisphere_lab,
                                   Brain_Cerebellum_lab,
                                   Brain_Cortex_lab,
                                   Brain_Frontal_Cortex_BA9_lab,
                                   Brain_Hippocampus_lab,
                                   Brain_Hypothalamus_lab,
                                   Brain_Nucleus_accumbens_basal_ganglia_lab,
                                   Brain_Putamen_basal_ganglia_lab,
                                   Brain_Spinal_cord_cervical_lab,
                                   Brain_Substantia_nigra_lab,
                                   Breast_Mammary_Tissue_lab,
                                   Colon_Sigmoid_lab,
                                   Colon_Transverse_lab,
                                   Esophagus_Gastroesophageal_Junction_lab,
                                   Esophagus_Mucosa_lab,
                                   Esophagus_Muscularis_lab,
                                   Heart_Atrial_Appendage_lab,
                                   Heart_Left_Ventricle_lab,
                                   Kidney_Cortex_lab,
                                   Liver_lab,
                                   Lung_lab,
                                   Minor_Salivary_Gland_lab,
                                   Muscle_Skeletal_lab,
                                   Nerve_Tibial_lab,
                                   Ovary_lab,
                                   Pancreas_lab,
                                   Pituitary_lab,
                                   Skin_Not_Sun_Exposed_Suprapubic_lab,
                                   Skin_Sun_Exposed_Lower_leg_lab,
                                   Small_Intestine_Terminal_Ileum_lab,
                                   Spleen_lab,
                                   Stomach_lab,
                                   Thyroid_lab,
                                   Uterus_lab,
                                   Vagina_lab,
                                   Whole_Blood_lab,
                                   placenta_lab),
                       type = c(gtex_lab, placenta_lab))

chr8_data$tissues = factor(chr8_data$tissues, levels = unique(chr8_data$tissues))

png("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_s8.png", width = 11, height = 13, units = "in", res = 300)
ggplot(chr8_data, aes(x=tissues, y=median, fill=type)) +
  geom_violin() +
  theme_bw() +
  labs(y='Median allele balance', title = 'chr8', x='Tissues') +
  theme(legend.title=element_blank(), legend.position="top", plot.title = element_text(hjust = 0.5, size=16), axis.title = element_text(size=16), legend.text = element_text(margin = margin(r = 30, unit = "pt"), size=14), axis.text = element_text(size = 14)) +
  coord_cartesian(xlim=c(0.5, 1)) +
  geom_hline(yintercept = 0.5, linetype=2, col="black") +
  geom_hline(yintercept = 0.8, linetype=2, col="black") +
  coord_flip() +
  scale_fill_manual(values=c("#00BFC4", "#C77CFF"))
dev.off()
