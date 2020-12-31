library(ggplot2)
# Plot stacked plot
data = read.csv("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/chrX_escaping_samples_prop_per_gene.tsv", sep = "\t")

data_sort = data[order(-data$prop), ]
data_sort_2 = data_sort[order(data_sort$status), ]

data_sort_2$gene <- as.character(data_sort_2$gene)
data_sort_2$gene <- factor(data_sort_2$gene, levels=unique(data_sort_2$gene))

png("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_s9.png", width = 13, height=6, units = "in", res=300)
ggplot(data_sort_2, aes(x = gene, y = prop)) + 
  geom_bar(stat = "identity", aes(fill=status)) +
  xlab("\nGene") +
  ylab("Proportion of samples") +
  theme_bw() +
  scale_fill_manual(values=c("navyblue", "gold"),
                    labels=c("Escape", "Inactivated")) +
  theme(axis.text.x=element_blank(), axis.text.y = element_text(size = 14), axis.title=element_text(size=16), legend.text=element_text(size=16), legend.title=element_text(size=18)) +
  geom_vline(xintercept = 22.5, size=1.5) +
  geom_vline(xintercept = 54.5, size=1.5)
dev.off()