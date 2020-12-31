library(ggplot2)
library(ggpubr)

# gtex tissue
# chrX
setwd('/scratch/tphung3/Placenta_XCI/gtex/02_analyze_gtex_counts/chrX_WASP_raw_counts_allele_balance/')
chrX_Adipose_Subcutaneous = read.csv('Adipose_Subcutaneous_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Adipose_Visceral_Omentum = read.csv('Adipose_Visceral_Omentum_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Adrenal_Gland = read.csv('Adrenal_Gland_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Artery_Aorta = read.csv('Artery_Aorta_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Artery_Coronary = read.csv('Artery_Coronary_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Artery_Tibial = read.csv('Artery_Tibial_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Amygdala = read.csv('Brain_Amygdala_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Anterior_cingulate_cortex_BA24 = read.csv('Brain_Anterior_cingulate_cortex_BA24_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Caudate_basal_ganglia = read.csv('Brain_Caudate_basal_ganglia_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Cerebellar_Hemisphere = read.csv('Brain_Cerebellar_Hemisphere_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Cerebellum = read.csv('Brain_Cerebellum_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Cortex = read.csv('Brain_Cortex_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Frontal_Cortex_BA9 = read.csv('Brain_Frontal_Cortex_BA9_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Hippocampus = read.csv('Brain_Hippocampus_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Hypothalamus = read.csv('Brain_Hypothalamus_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Nucleus_accumbens_basal_ganglia = read.csv('Brain_Nucleus_accumbens_basal_ganglia_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Putamen_basal_ganglia = read.csv('Brain_Putamen_basal_ganglia_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Spinal_cord_cervical = read.csv('Brain_Spinal_cord_cervical_c-1_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Brain_Substantia_nigra = read.csv('Brain_Substantia_nigra_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Breast_Mammary_Tissue = read.csv('Breast_Mammary_Tissue_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Colon_Sigmoid = read.csv('Colon_Sigmoid_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Colon_Transverse = read.csv('Colon_Transverse_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Esophagus_Gastroesophageal_Junction = read.csv('Esophagus_Gastroesophageal_Junction_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Esophagus_Mucosa = read.csv('Esophagus_Mucosa_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Esophagus_Muscularis = read.csv('Esophagus_Muscularis_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Heart_Atrial_Appendage = read.csv('Heart_Atrial_Appendage_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Heart_Left_Ventricle = read.csv('Heart_Left_Ventricle_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Kidney_Cortex = read.csv('Kidney_Cortex_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Liver = read.csv('Liver_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Lung = read.csv('Lung_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Minor_Salivary_Gland = read.csv('Minor_Salivary_Gland_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Muscle_Skeletal = read.csv('Muscle_Skeletal_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Nerve_Tibial = read.csv('Nerve_Tibial_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Ovary = read.csv('Ovary_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Pancreas = read.csv('Pancreas_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Pituitary = read.csv('Pituitary_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Skin_Not_Sun_Exposed_Suprapubic = read.csv('Skin_Not_Sun_Exposed_Suprapubic_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Skin_Sun_Exposed_Lower_leg = read.csv('Skin_Sun_Exposed_Lower_leg_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Small_Intestine_Terminal_Ileum = read.csv('Small_Intestine_Terminal_Ileum_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Spleen = read.csv('Spleen_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Stomach = read.csv('Stomach_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Thyroid = read.csv('Thyroid_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Uterus = read.csv('Uterus_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Vagina = read.csv('Vagina_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_Whole_Blood = read.csv('Whole_Blood_chrX_nonpars_median_allele_balance.txt', sep = '\t')


# Placenta tissues
chrX_placenta = read.csv('/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chrX/placenta_chrX_nonpars_median_allele_balance.txt', sep = '\t')
chrX_placenta_rm = subset(chrX_placenta, chrX_placenta$dna_id != "OBG0175")

Adipose_Subcutaneous_lab = rep("Adipose Subcutaneous", nrow(chrX_Adipose_Subcutaneous))
Adipose_Visceral_Omentum_lab = rep("Adipose Visceral Omentum", nrow(chrX_Adipose_Visceral_Omentum))
Adrenal_Gland_lab = rep("Adrenal Gland", nrow(chrX_Adrenal_Gland))
Artery_Aorta_lab = rep("Artery Aorta", nrow(chrX_Artery_Aorta))
Artery_Coronary_lab = rep("Artery Coronary", nrow(chrX_Artery_Coronary))
Artery_Tibial_lab = rep("Artery Tibial", nrow(chrX_Artery_Tibial))
Brain_Amygdala_lab = rep("Brain Amygdala", nrow(chrX_Brain_Amygdala))
Brain_Anterior_cingulate_cortex_BA24_lab = rep("Brain Anterior cingulate cortex BA24", nrow(chrX_Brain_Anterior_cingulate_cortex_BA24))
Brain_Caudate_basal_ganglia_lab = rep("Brain Caudate basal ganglia", nrow(chrX_Brain_Caudate_basal_ganglia))
Brain_Cerebellar_Hemisphere_lab = rep("Brain Cerebellar Hemisphere", nrow(chrX_Brain_Cerebellar_Hemisphere))
Brain_Cerebellum_lab = rep("Brain Cerebellum", nrow(chrX_Brain_Cerebellum))
Brain_Cortex_lab = rep("Brain Cortex", nrow(chrX_Brain_Cortex))
Brain_Frontal_Cortex_BA9_lab = rep("Brain Frontal Cortex BA9", nrow(chrX_Brain_Frontal_Cortex_BA9))
Brain_Hippocampus_lab = rep("Brain Hippocampus", nrow(chrX_Brain_Hippocampus))
Brain_Hypothalamus_lab = rep("Brain Hypothalamus", nrow(chrX_Brain_Hypothalamus))
Brain_Nucleus_accumbens_basal_ganglia_lab = rep("Brain Nucleus accumbens basal ganglia", nrow(chrX_Brain_Nucleus_accumbens_basal_ganglia))
Brain_Putamen_basal_ganglia_lab = rep("Brain Putamen basal ganglia", nrow(chrX_Brain_Putamen_basal_ganglia))
Brain_Spinal_cord_cervical_lab = rep("Brain Spinal cord cervical", nrow(chrX_Brain_Spinal_cord_cervical))
Brain_Substantia_nigra_lab = rep("Brain Substantia nigra", nrow(chrX_Brain_Substantia_nigra))
Breast_Mammary_Tissue_lab = rep("Breast Mammary Tissue", nrow(chrX_Breast_Mammary_Tissue))
Colon_Sigmoid_lab = rep("Colon Sigmoid", nrow(chrX_Colon_Sigmoid))
Colon_Transverse_lab = rep("Colon Transverse", nrow(chrX_Colon_Transverse))
Esophagus_Gastroesophageal_Junction_lab = rep("Esophagus Gastroesophageal Junction", nrow(chrX_Esophagus_Gastroesophageal_Junction))
Esophagus_Mucosa_lab = rep("Esophagus Mucosa", nrow(chrX_Esophagus_Mucosa))
Esophagus_Muscularis_lab = rep("Esophagus Muscularis", nrow(chrX_Esophagus_Muscularis))
Heart_Atrial_Appendage_lab = rep("Heart Atrial Appendage", nrow(chrX_Heart_Atrial_Appendage))
Heart_Left_Ventricle_lab = rep("Heart Left Ventricle", nrow(chrX_Heart_Left_Ventricle))
Kidney_Cortex_lab = rep("Kidney Cortex", nrow(chrX_Kidney_Cortex))
Liver_lab = rep("Liver", nrow(chrX_Liver))
Lung_lab = rep("Lung", nrow(chrX_Lung))
Minor_Salivary_Gland_lab = rep("Minor Salivary Gland", nrow(chrX_Minor_Salivary_Gland))
Muscle_Skeletal_lab = rep("Muscle Skeletal", nrow(chrX_Muscle_Skeletal))
Nerve_Tibial_lab = rep("Nerve Tibial", nrow(chrX_Nerve_Tibial))
Ovary_lab = rep("Ovary", nrow(chrX_Ovary))
Pancreas_lab = rep("Pancreas", nrow(chrX_Pancreas))
Pituitary_lab = rep("Pituitary", nrow(chrX_Pituitary))
Skin_Not_Sun_Exposed_Suprapubic_lab = rep("Skin Not Sun Exposed Suprapubic", nrow(chrX_Skin_Not_Sun_Exposed_Suprapubic))
Skin_Sun_Exposed_Lower_leg_lab = rep("Skin Sun Exposed Lower leg", nrow(chrX_Skin_Sun_Exposed_Lower_leg))
Small_Intestine_Terminal_Ileum_lab = rep("Small Intestine Terminal Ileum", nrow(chrX_Small_Intestine_Terminal_Ileum))
Spleen_lab = rep("Spleen", nrow(chrX_Spleen))
Stomach_lab = rep("Stomach", nrow(chrX_Stomach))
Thyroid_lab = rep("Thyroid", nrow(chrX_Thyroid))
Uterus_lab = rep("Uterus", nrow(chrX_Uterus))
Vagina_lab = rep("Vagina", nrow(chrX_Vagina))
Whole_Blood_lab = rep("Whole Blood", nrow(chrX_Whole_Blood))

placenta_lab = rep('Placenta', 58)

n_gtex = nrow(chrX_Thyroid)+ nrow(chrX_Uterus)+ nrow(chrX_Pancreas)+ nrow(chrX_Artery_Aorta)+ nrow(chrX_Spleen)+ nrow(chrX_Small_Intestine_Terminal_Ileum)+ nrow(chrX_Skin_Not_Sun_Exposed_Suprapubic)+ nrow(chrX_Pituitary)+ nrow(chrX_Esophagus_Mucosa)+ nrow(chrX_Brain_Cortex)+ nrow(chrX_Artery_Tibial)+ nrow(chrX_Artery_Coronary)+ nrow(chrX_Kidney_Cortex)+ nrow(chrX_Nerve_Tibial)+ nrow(chrX_Brain_Cerebellum)+ nrow(chrX_Lung)+ nrow(chrX_Adrenal_Gland)+ nrow(chrX_Heart_Left_Ventricle)+ nrow(chrX_Brain_Caudate_basal_ganglia)+ nrow(chrX_Esophagus_Muscularis)+ nrow(chrX_Esophagus_Gastroesophageal_Junction)+ nrow(chrX_Brain_Cerebellar_Hemisphere)+ nrow(chrX_Brain_Spinal_cord_cervical)+ nrow(chrX_Breast_Mammary_Tissue)+ nrow(chrX_Brain_Frontal_Cortex_BA9)+ nrow(chrX_Brain_Hypothalamus)+ nrow(chrX_Minor_Salivary_Gland)+ nrow(chrX_Whole_Blood)+ nrow(chrX_Brain_Substantia_nigra)+ nrow(chrX_Liver)+ nrow(chrX_Colon_Transverse)+ nrow(chrX_Heart_Atrial_Appendage)+ nrow(chrX_Brain_Amygdala)+ nrow(chrX_Stomach)+ nrow(chrX_Ovary)+ nrow(chrX_Adipose_Visceral_Omentum)+ nrow(chrX_Skin_Sun_Exposed_Lower_leg)+ nrow(chrX_Colon_Sigmoid)+ nrow(chrX_Brain_Anterior_cingulate_cortex_BA24)+ nrow(chrX_Adipose_Subcutaneous)+ nrow(chrX_Brain_Hippocampus)+ nrow(chrX_Brain_Nucleus_accumbens_basal_ganglia)+ nrow(chrX_Brain_Putamen_basal_ganglia)+ nrow(chrX_Vagina)+ nrow(chrX_Muscle_Skeletal)

gtex_lab = rep('GTEX', n_gtex)

chrX_data = data.frame(median = c(chrX_Adipose_Subcutaneous$median_allele_balance_subset,
                                  chrX_Adipose_Visceral_Omentum$median_allele_balance_subset,
                                  chrX_Adrenal_Gland$median_allele_balance_subset,
                                  chrX_Artery_Aorta$median_allele_balance_subset,
                                  chrX_Artery_Coronary$median_allele_balance_subset,
                                  chrX_Artery_Tibial$median_allele_balance_subset,
                                  chrX_Brain_Amygdala$median_allele_balance_subset,
                                  chrX_Brain_Anterior_cingulate_cortex_BA24$median_allele_balance_subset,
                                  chrX_Brain_Caudate_basal_ganglia$median_allele_balance_subset,
                                  chrX_Brain_Cerebellar_Hemisphere$median_allele_balance_subset,
                                  chrX_Brain_Cerebellum$median_allele_balance_subset,
                                  chrX_Brain_Cortex$median_allele_balance_subset,
                                  chrX_Brain_Frontal_Cortex_BA9$median_allele_balance_subset,
                                  chrX_Brain_Hippocampus$median_allele_balance_subset,
                                  chrX_Brain_Hypothalamus$median_allele_balance_subset,
                                  chrX_Brain_Nucleus_accumbens_basal_ganglia$median_allele_balance_subset,
                                  chrX_Brain_Putamen_basal_ganglia$median_allele_balance_subset,
                                  chrX_Brain_Spinal_cord_cervical$median_allele_balance_subset,
                                  chrX_Brain_Substantia_nigra$median_allele_balance_subset,
                                  chrX_Breast_Mammary_Tissue$median_allele_balance_subset,
                                  chrX_Colon_Sigmoid$median_allele_balance_subset,
                                  chrX_Colon_Transverse$median_allele_balance_subset,
                                  chrX_Esophagus_Gastroesophageal_Junction$median_allele_balance_subset,
                                  chrX_Esophagus_Mucosa$median_allele_balance_subset,
                                  chrX_Esophagus_Muscularis$median_allele_balance_subset,
                                  chrX_Heart_Atrial_Appendage$median_allele_balance_subset,
                                  chrX_Heart_Left_Ventricle$median_allele_balance_subset,
                                  chrX_Kidney_Cortex$median_allele_balance_subset,
                                  chrX_Liver$median_allele_balance_subset,
                                  chrX_Lung$median_allele_balance_subset,
                                  chrX_Minor_Salivary_Gland$median_allele_balance_subset,
                                  chrX_Muscle_Skeletal$median_allele_balance_subset,
                                  chrX_Nerve_Tibial$median_allele_balance_subset,
                                  chrX_Ovary$median_allele_balance_subset,
                                  chrX_Pancreas$median_allele_balance_subset,
                                  chrX_Pituitary$median_allele_balance_subset,
                                  chrX_Skin_Not_Sun_Exposed_Suprapubic$median_allele_balance_subset,
                                  chrX_Skin_Sun_Exposed_Lower_leg$median_allele_balance_subset,
                                  chrX_Small_Intestine_Terminal_Ileum$median_allele_balance_subset,
                                  chrX_Spleen$median_allele_balance_subset,
                                  chrX_Stomach$median_allele_balance_subset,
                                  chrX_Thyroid$median_allele_balance_subset,
                                  chrX_Uterus$median_allele_balance_subset,
                                  chrX_Vagina$median_allele_balance_subset,
                                  chrX_Whole_Blood$median_allele_balance_subset,
                                  chrX_placenta_rm$median_allele_balance_subset),
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

chrX_data$tissues = factor(chrX_data$tissues, levels = unique(chrX_data$tissues))

png("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_3.png", width = 11, height = 13, units = "in", res = 300)
ggplot(chrX_data, aes(x=tissues, y=median, fill=type)) +
  geom_violin() +
  theme_bw() +
  labs(y='Median allele balance', title = 'chrX', x='Tissues') +
  theme(legend.title=element_blank(), legend.position="top", plot.title = element_text(hjust = 0.5, size=16), axis.title = element_text(size=16), legend.text = element_text(margin = margin(r = 30, unit = "pt"), size=14), axis.text = element_text(size = 14)) +
  coord_cartesian(xlim=c(0.5, 1)) +
  geom_hline(yintercept = 0.5, linetype=2, col="black") +
  geom_hline(yintercept = 0.8, linetype=2, col="black") +
  scale_fill_manual(values=c("#00BFC4", "#C77CFF")) +
  coord_flip() 
dev.off()