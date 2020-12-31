library(ggplot2)

setwd("/scratch/tphung3/Placenta_XCI/gene_level/female_male_log2ratio/")

inactivated = read.table("inactivated_genes_log2ratio.tsv")
escape = read.table("escape_genes_log2ratio.tsv")
variable = read.table("variable_genes_log2ratio.tsv")

inactivated_labels = rep("Inactivated", nrow(inactivated))
escape_labels = rep("Escape", nrow(escape))
variable_labels = rep("Variable", nrow(variable))

data = data.frame(log2ratio = c(inactivated$V3, escape$V3, variable$V3), labs = c(inactivated_labels, escape_labels, variable_labels))

data$labs <- as.character(data$labs)
data$labs <- factor(data$labs, levels=unique(data$labs))

png('/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_sx_log2ratio.png', width = 9, height = 6, units = "in", res = 300)
ggplot(data, aes(x=labs, y=log2ratio, color=labs)) +
  geom_violin(position = position_dodge(width = 0.9)) + 
  geom_boxplot(width = 0.1, outlier.shape = NA) + 
  geom_jitter() +
  theme_bw() + 
  labs(x="Inactivation status", y="Female/Male Log2Ratio") +
  scale_color_manual(name = "", labels = c("Inactivated", "Escape", "Variable"), values = c("gold", "navyblue", "skyblue")) +
  theme(axis.text = element_text(color = "black", size=14), axis.title=element_text(size = 16), legend.position="none") 
dev.off()