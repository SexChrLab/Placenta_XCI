library(ggplot2)

setwd('/scratch/tphung3/PlacentaSexDiff/E_escape_genes/gtex_counts/genes_allele_balance/chrX/')

args = commandArgs(trailingOnly=TRUE)

gene = args[1]

gtex = read.csv(paste(gene, '_allele_balance_for_gtex_samples.tsv', sep=''), sep = '\t')
placenta = read.csv(paste(gene, '_allele_balance_for_placenta_samples.tsv', sep=''), sep = '\t')
decidua_females = read.csv(paste(gene, '_allele_balance_for_decidua_females_samples.tsv', sep=''), sep = '\t')
decidua_males = read.csv(paste(gene, '_allele_balance_for_decidua_males_samples.tsv', sep=''), sep = '\t')

gtex_label = rep('GTEX', nrow(gtex))
placenta_label = rep('Placenta', nrow(placenta))
decidua_females_label = rep('Decidua Females', nrow(decidua_females))
decidua_males_label = rep('Decidua Males', nrow(decidua_males))

data = data.frame(allele_balance = c(gtex$allele_balance, placenta$allele_balance, decidua_females$allele_balance, decidua_males$allele_balance), labels = c(gtex_label, placenta_label, decidua_females_label, decidua_males_label))

data$labels = factor(data$labels, levels = unique(data$labels))

png(paste('/scratch/tphung3/PlacentaSexDiff/E_escape_genes/gtex_counts/genes_allele_balance/chrX/plots/', gene, '_allele_balance_compare_gtex_placenta_decidua_violin.png', sep=''), width = 11, height = 6, units = "in", res = 300)
ggplot(data, aes(x=data$labels, y=data$allele_balance)) +
  geom_violin() +
  theme_bw() +
  labs(x="Tissue", y="Allele balance", title = gene) +
  theme(plot.title = element_text(hjust = 0.5, size=16), axis.text = element_text(size = 14), axis.title = element_text(size = 14)) +
  coord_cartesian(ylim=c(0.5, 1)) +
  geom_hline(yintercept = 0.8, col="darkgray", linetype=2) +
  geom_hline(yintercept = 1, col="darkgray", linetype=2) +
  geom_jitter(shape=16, position=position_jitter(0.2))
dev.off()
