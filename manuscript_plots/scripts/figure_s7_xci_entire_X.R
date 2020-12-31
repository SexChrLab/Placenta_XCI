# ---------------------------------------
# For Figure SX
# Plot for the entire X chromosome phased
# ---------------------------------------
library(ggplot2)

# Set up directory paths
allele_balance_dir = "/scratch/tphung3/PlacentaSexDiff/A_placenta/08_phasing/phased_allele_balance/"
save_plot_dir = "/scratch/tphung3/PlacentaSexDiff/C_manuscript/manuscript_plots/"

# Write a function for plotting
plot_func = function(data, sample_id) {
  ggplot() +
    geom_rect(aes(xmin=10001, xmax=2781479, ymin=0, ymax=Inf), fill="gray", color="black") +
    geom_rect(aes(xmin=73820892, xmax=73851867, ymin=0, ymax=Inf), fill="gray", color="black") +
    geom_rect(aes(xmin=155701383, xmax=156030895, ymin=0, ymax=Inf), fill="gray", color="black") +
    geom_point(data, mapping = aes(x=position, y=phased_x, colour="SiteA", shape="SiteA", size="SiteA", fill="SiteA"), stroke=1.5) +
    geom_point(data, mapping = aes(x=position, y=phased_y, colour="SiteB", shape="SiteB", size="SiteB", fill="SiteB"), stroke=1.5) +
    coord_cartesian(xlim = c(1, 156040895), ylim=c(0,1))+
    scale_x_continuous(breaks=c(0, 20000000, 40000000, 60000000, 80000000, 100000000, 120000000, 140000000),
                       labels=c("0", "20", "40", "60", "80", "100", "120", "140")) +
    theme_bw() +
    labs(x="chrX position (Mb)", y="Allele balance", title=sample_id) +
    theme(plot.title = element_text(hjust = 0.5, size=22), axis.text.x = element_text(size=18, colour="black"), axis.text.y = element_text(size=18, colour = "black"), axis.title= element_text(size=20), panel.background = element_blank(), axis.line=element_line(colour = "black"), axis.ticks = element_line(colour = "black", size=1.5), axis.ticks.length = unit(0.3, 'cm'), legend.position="top", legend.title = element_blank(), legend.text = element_text(size=16)) +
    geom_hline(yintercept = 0.8, color="darkgray", size=0.75, linetype=2) +
    geom_hline(yintercept = 0.2, color="darkgray", size=0.75, linetype=2) +
    scale_colour_manual(name="Placenta", values=c(SiteA="black", SiteB="red")) +
    scale_shape_manual(name="Placenta", values=c(SiteA=21, SiteB=24)) +
    scale_size_manual(name="Placenta", values=c(SiteA=2, SiteB=1.2)) +
    scale_fill_manual(name="Placenta", values=c(SiteA=NA, SiteB="red"))
}

# Create list of data
samples <- list(
  "OBG0030"="OBG-0030-PLAC_chrX_phased_allele_balance_data.tsv",
  "OBG0201"="OBG-0201-PLAC_chrX_phased_allele_balance_data.tsv",
  "OBG0051"="OBG-0051-PLAC_chrX_phased_allele_balance_data.tsv"
)


# Plot and save
for (sample in names(samples)) {
  placenta_data = read.csv(paste(allele_balance_dir, samples[[sample]], sep=""))
  p = plot_func(data = placenta_data, sample_id = sample)
  ggsave(paste(save_plot_dir, "figure_sx_xci_entire_X_", sample, ".png"), plot = p, width = 11,height=4)
}