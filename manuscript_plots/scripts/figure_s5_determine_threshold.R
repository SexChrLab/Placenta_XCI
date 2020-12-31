library(ggplot2)
library(ggpubr)

args = commandArgs(trailingOnly=TRUE)
dna = args[1]
rna_1 = args[2]
rna_2 = args[3]

# site A
siteA_chrX = read.csv(paste("/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chrX/", dna, "_", rna_1, "_chrX_allele_balance.tsv", sep=""), sep="\t")
siteA_chr8 = read.csv(paste("/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chr8/", dna, "_", rna_1, "_chr8_allele_balance.tsv", sep=""), sep="\t")

siteA_chrX_filtered = subset(siteA_chrX, siteA_chrX$total_count>10)
siteA_chr8_filtered = subset(siteA_chr8, siteA_chr8$total_count>10)

siteA_df = data.frame(
  Chromosome=c(rep("chrX", nrow(siteA_chrX_filtered)), rep("chr8", nrow(siteA_chr8_filtered))),
  allele_balance=c(siteA_chrX_filtered$allele_balance, siteA_chr8_filtered$allele_balance)
)

p1 = ggplot(siteA_df, aes(x=allele_balance, color=Chromosome, fill=Chromosome)) +
  geom_histogram(alpha=0.5, position="identity") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme_bw() +
  labs(x="Allele balance", y="Count", title = "Site A") +
  theme(axis.title = element_text(size=12),
        plot.title = element_text(hjust = 0.5, size=16),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        legend.spacing.x = unit(0.2, 'cm')) +
  geom_vline(xintercept = 0.8, linetype=2)

# site B
siteB_chrX = read.csv(paste("/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chrX/", dna, "_", rna_2, "_chrX_allele_balance.tsv", sep=""), sep="\t")
siteB_chr8 = read.csv(paste("/scratch/tphung3/Placenta_XCI/placenta/03_analyze_ase/results/chr8/", dna, "_", rna_2, "_chr8_allele_balance.tsv", sep=""), sep="\t")

siteB_chrX_filtered = subset(siteB_chrX, siteB_chrX$total_count>10)
siteB_chr8_filtered = subset(siteB_chr8, siteB_chr8$total_count>10)

siteB_df = data.frame(
  Chromosome=c(rep("chrX", nrow(siteB_chrX_filtered)), rep("chr8", nrow(siteB_chr8_filtered))),
  allele_balance=c(siteB_chrX_filtered$allele_balance, siteB_chr8_filtered$allele_balance)
)

p2 = ggplot(siteB_df, aes(x=allele_balance, color=Chromosome, fill=Chromosome)) +
  geom_histogram(alpha=0.5, position="identity") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  theme_bw() +
  labs(x="Allele balance", y="Count", title = "Site B") +
  theme(axis.title = element_text(size=12),
        plot.title = element_text(hjust = 0.5, size=16),
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.spacing.x = unit(0.2, 'cm')) +
  geom_vline(xintercept = 0.8, linetype=2)

png(paste("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_s5/", dna, ".png", sep=""), width = 10, height = 4, units = "in", res = 300)
ggarrange(p1, p2, ncol = 2, common.legend = T)
dev.off()