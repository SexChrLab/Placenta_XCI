library(ggplot2)
library(ggpubr)
setwd("/scratch/tphung3/Placenta_XCI/placenta/06_nonPAR_males/asereadcounter/")
dir.rl = getwd()

samples <- list(
  "OBG0112"=c('OBG0112/OBG0112_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'OBG0112/OBG0112_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv'),
  "OBG0116"=c('OBG0116/OBG0116_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'OBG0116/OBG0116_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv'),
  "OBG0117"=c('OBG0117/OBG0117_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'OBG0117/OBG0117_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv'),
  "OBG0118"=c('OBG0118/OBG0118_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'OBG0118/OBG0118_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv'),
  "OBG0122"=c('OBG0122/OBG0122_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'OBG0122/OBG0122_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv'),
  "OBG0123"=c('OBG0123/OBG0123_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'OBG0123/OBG0123_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv'),
  "OBG0126"=c('OBG0126/OBG0126_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'OBG0126/OBG0126_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv'),
  "OBG0130"=c('OBG0130/OBG0130_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'OBG0130/OBG0130_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv'),
  "OBG0132"=c('OBG0132/OBG0132_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'OBG0132/OBG0132_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv'),
  "OBG0158"=c('OBG0158/OBG0158_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'OBG0158/OBG0158_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv'),
  "YPOPS0006"=c('YPOPS0006/YPOPS0006_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'YPOPS0006/YPOPS0006_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv'),
  "OBG0053"=c('OBG0053/OBG0053_placenta_1_nonPAR_diploid_filter_for_biallelic_hets.tsv', 'OBG0053/OBG0053_placenta_2_nonPAR_diploid_filter_for_biallelic_hets.tsv')
)

for (sample in names(samples)) {
  # placenta_1
  placenta_1_chrX = read.csv(samples[[sample]][1], sep = '\t')
  colnames(placenta_1_chrX) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
  
  placenta_1_chrX$refRatio = placenta_1_chrX$refCount/placenta_1_chrX$totalCount
  placenta_1_chrX$altRatio = placenta_1_chrX$altCount/placenta_1_chrX$totalCount
  
  
  placenta_1_chrX_ratio = c()
  for (i in 1:nrow(placenta_1_chrX)) {
    row = placenta_1_chrX[i,]
    if (row$refRatio > row$altRatio) {
      placenta_1_chrX_ratio = c(placenta_1_chrX_ratio, row$refRatio)
    }
    else (placenta_1_chrX_ratio = c(placenta_1_chrX_ratio, row$altRatio))
  }
  placenta_1_chrX$ratio = placenta_1_chrX_ratio
  
  p1 = ggplot(placenta_1_chrX, aes(ratio)) +
    geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x="Allele balance", y="Count", title="Site A") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=18), axis.text.x = element_text(size=16, colour="black"), axis.text.y = element_text(size=16, colour = "black"), axis.title= element_text(size=16), panel.background = element_blank(), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), legend.title = element_blank(), legend.text = element_text(size=16), legend.position = "top") +
    geom_vline(xintercept = 0.8, color="darkgray", size=0.75, linetype=2) +
    geom_vline(xintercept = 0.2, color="darkgray", size=0.75, linetype=2)
  
  # placenta_2
  placenta_2_chrX = read.csv(samples[[sample]][2], sep = '\t')
  colnames(placenta_2_chrX) = c('contig', 'position', 'variantID', 'refAllele', 'altAllele', 'refCount', 'altCount', 'totalCount', 'lowMAPQDepth', 'lowBaseQDepth', 'rawDepth', 'otherBases', 'improperPairs')
  
  placenta_2_chrX$refRatio = placenta_2_chrX$refCount/placenta_2_chrX$totalCount
  placenta_2_chrX$altRatio = placenta_2_chrX$altCount/placenta_2_chrX$totalCount
  
  
  placenta_2_chrX_ratio = c()
  for (i in 1:nrow(placenta_2_chrX)) {
    row = placenta_2_chrX[i,]
    if (row$refRatio > row$altRatio) {
      placenta_2_chrX_ratio = c(placenta_2_chrX_ratio, row$refRatio)
    }
    else (placenta_2_chrX_ratio = c(placenta_2_chrX_ratio, row$altRatio))
  }
  placenta_2_chrX$ratio = placenta_2_chrX_ratio
  
  p2 = ggplot(placenta_2_chrX, aes(ratio)) +
    geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') +
    coord_cartesian(xlim = c(0, 1)) +
    labs(x="Allele balance", y="Count", title="Site B") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size=18), axis.text.x = element_text(size=16, colour="black"), axis.text.y = element_text(size=16, colour = "black"), axis.title= element_text(size=16), panel.background = element_blank(), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), legend.title = element_blank(), legend.text = element_text(size=16), legend.position = "top") +
    geom_vline(xintercept = 0.8, color="darkgray", size=0.75, linetype=2) +
    geom_vline(xintercept = 0.2, color="darkgray", size=0.75, linetype=2)
  gg = ggarrange(p1, p2, ncol = 2, nrow=1, common.legend = T)
  ggsave(file.path("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_s6/",sprintf("figure_s6_%s.png",sample)),plot = gg, width = 10,height=4, dpi=300)
}