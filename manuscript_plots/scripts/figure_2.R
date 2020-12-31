library(ggplot2)
library(ggpubr)

# Set working directory. Change this here if you need to rerun.
setwd("/scratch/tphung3/Placenta_XCI/placenta/04_phasing/phased_allele_balance/") 

# chrX
chrX = read.table("all_placenta_chrX_phased_allele_balance.tsv")
# Generate a dataframe for plotting
# Note: In the samples columns, need to make sure that the orders are correct. I manually make this so that the naming convention of the samples is the same. 
chrX_df_plot = data.frame(vals = c(chrX[,2], chrX[,3]), Sites = c("Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B"), samples = c("OBG0170", "OBG0174", "OBG0133", "OBG0120", "OBG0068", "OBG0342", "OBG0024", "OBG0030", "OBG0188", "OBG0289", "OBG0338", "OBG0205", "OBG0026", "OBG0066", "OBG0111", "OBG0051", "OBG0039", "OBG0156", "OBG0180", "OBG0022", "OBG0175", "OBG0044", "OBG0121", "OBG0178", "OBG0050", "OBG0115", "OBG0166", "OBG0028", "OBG0201", "OBG0138", "OBG0170", "OBG0174", "OBG0133", "OBG0120", "OBG0068", "OBG0342", "OBG0024", "OBG0030", "OBG0188", "OBG0289", "OBG0338", "OBG0205", "OBG0026", "OBG0066", "OBG0111", "OBG0051", "OBG0039", "OBG0156", "OBG0180", "OBG0022", "OBG0175", "OBG0044", "OBG0121", "OBG0178", "OBG0050", "OBG0115", "OBG0166", "OBG0028", "OBG0201", "OBG0138"))

chrX_df_plot$samples <- as.character(chrX_df_plot$samples)
chrX_df_plot$samples <- factor(chrX_df_plot$samples, levels=unique(chrX_df_plot$samples))

# Remove the sample OBG0175
chrX_df_plot_rm = subset(chrX_df_plot, chrX_df_plot$samples != "OBG0175")

# Change the output path here when rerun
png("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_2_phased_allele_balance_chrX.png", width = 11, height = 5, units = "in", res = 1200)
ggplot(chrX_df_plot_rm, aes(x=samples, y=vals, color=Sites, shape=Sites, size=Sites, fill=Sites)) +
  geom_point(stroke=2) +
  scale_colour_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(21, 24)) +
  scale_size_manual(values=c(7, 5)) +
  scale_fill_manual(values=c(NA, "red")) +
  labs(x="", y="Allele balance") +
  theme(axis.text.x = element_text(size=16, angle = 60, hjust=1), axis.text.y = element_text(size=16), axis.title= element_text(size=18), panel.background = element_blank(), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), legend.text = element_text(size=16), legend.title = element_text(size = 18), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.title = element_text(hjust = 0.5, size = 24)) +
  geom_hline(yintercept = 0.8, color="darkgray", size=0.75, linetype=2) +
  geom_hline(yintercept = 0.2, color="darkgray", size=0.75, linetype=2)
dev.off()

# chr8
chr8 = read.table("all_placenta_chr8_phased_allele_balance.tsv")
chr8_modified = chr8[match(chrX[,1], chr8[,1]),]

chr8_df_plot = data.frame(vals = c(chr8_modified[,2], chr8_modified[,3]), Sites = c("Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site A", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B", "Site B"), samples = c("OBG0170", "OBG0174", "OBG0133", "OBG0120", "OBG0068", "OBG0342", "OBG0024", "OBG0030", "OBG0188", "OBG0289", "OBG0338", "OBG0205", "OBG0026", "OBG0066", "OBG0111", "OBG0051", "OBG0039", "OBG0156", "OBG0180", "OBG0022", "OBG0175", "OBG0044", "OBG0121", "OBG0178", "OBG0050", "OBG0115", "OBG0166", "OBG0028", "OBG0201", "OBG0138", "OBG0170", "OBG0174", "OBG0133", "OBG0120", "OBG0068", "OBG0342", "OBG0024", "OBG0030", "OBG0188", "OBG0289", "OBG0338", "OBG0205", "OBG0026", "OBG0066", "OBG0111", "OBG0051", "OBG0039", "OBG0156", "OBG0180", "OBG0022", "OBG0175", "OBG0044", "OBG0121", "OBG0178", "OBG0050", "OBG0115", "OBG0166", "OBG0028", "OBG0201", "OBG0138"))

chr8_df_plot$samples <- as.character(chr8_df_plot$samples)
chr8_df_plot$samples <- factor(chr8_df_plot$samples, levels=unique(chr8_df_plot$samples))

# Remove the sample OBG0175
chr8_df_plot_rm = subset(chr8_df_plot, chr8_df_plot$samples != "OBG0175")

png("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_2_phased_allele_balance_chr8.png", width = 11, height = 5, units = "in", res = 300)
ggplot(chr8_df_plot, aes(x=samples, y=vals, color=Sites, shape=Sites, size=Sites, fill=Sites)) +
  geom_point(stroke=2) +
  scale_colour_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(21, 24)) +
  scale_size_manual(values=c(7, 5)) +
  scale_fill_manual(values=c(NA, "red")) +
  labs(x="", y="Allele balance") +
  theme(axis.text.x = element_text(size=16, angle = 60, hjust=1), axis.text.y = element_text(size=16), axis.title= element_text(size=18), panel.background = element_blank(), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), legend.text = element_text(size=16), legend.title = element_text(size = 18), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.title = element_text(hjust = 0.5, size = 24)) +
  geom_hline(yintercept = 0.8, color="darkgray", size=0.75, linetype=2) +
  geom_hline(yintercept = 0.2, color="darkgray", size=0.75, linetype=2) +
  coord_cartesian(ylim=c(0, 1))
dev.off()