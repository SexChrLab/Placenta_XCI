library(ggplot2)
# chrX
chrX = read.table("/scratch/tphung3/Placenta_XCI/heart/phased_allele_balance/all_heart_chrX_phased_allele_balance.tsv")

chrX_df_plot = data.frame(vals = c(chrX[,3], chrX[,2]), Sites = c(rep("Heart Left Ventricle", 85), rep("Heart Atrial Appendage", 85)), samples = chrX[,1])

chrX_df_plot$samples <- as.character(chrX_df_plot$samples)
chrX_df_plot$samples <- factor(chrX_df_plot$samples, levels=unique(chrX_df_plot$samples))

png("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_2C.png", width = 11, height = 5, units = "in", res = 1200)
ggplot(chrX_df_plot, aes(x=samples, y=vals, color=Sites, shape=Sites, size=Sites, fill=Sites)) +
  geom_point(stroke=2) +
  scale_colour_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(21, 24)) +
  scale_size_manual(values=c(7, 5)) +
  scale_fill_manual(values=c(NA, "red")) +
  labs(x="Subjects", y="Allele balance") +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=16), axis.title= element_text(size=18), panel.background = element_blank(), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), legend.text = element_text(size=16), legend.title = element_text(size = 18), plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.title = element_text(hjust = 0.5, size = 24), legend.position = "top") +
  geom_hline(yintercept = 0.8, color="darkgray", size=0.75, linetype=2) +
  geom_hline(yintercept = 0.2, color="darkgray", size=0.75, linetype=2)
dev.off()