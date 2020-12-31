# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate")

# SNPrelate PCA
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(ggpubr)

setwd("/scratch/tphung3/Placenta_XCI/placenta/05_pca/")

pop_code = c(rep('Batch 1', 18), rep('Batch 2', 12))

### chrX all variants (before any filtering)
vcf.fn = "chrX.gatk.called.vqsr.sv.placentas.vcf"
snpgdsVCF2GDS(vcf.fn, "chrX_all_variants.gds", method="biallelic.only")
snpgdsSummary("chrX_all_variants.gds")
genofile <- snpgdsOpen("chrX_all_variants.gds")
pca = snpgdsPCA(genofile, autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# Make a data frame
sample.id = pca$sample.id
chrX_df <- data.frame(sample.id = pca$sample.id,
                      pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,2],    # the second eigenvector
                      EV3 = pca$eigenvect[,3],
                      EV4 = pca$eigenvect[,4],
                      EV5 = pca$eigenvect[,5],
                      stringsAsFactors = FALSE)


p1 = ggplot(chrX_df, aes(x=EV1, y=EV2, fill=pop, shape=pop)) + 
  geom_point(size=3, colour="black") + 
  theme_bw()+
  labs(title="X chromosome", x="PC1", y="PC2") +
  theme(axis.text.x = element_text(size=14, colour="black"), axis.text.y = element_text(size=14, colour = "black"), axis.title.y = element_text(size=18, angle=0, vjust=0.5), axis.title.x = element_text(size=18), axis.line=element_line(colour = "black"),plot.title = element_text(hjust = 0.5, size=18), legend.title = element_blank(), panel.background = element_blank(), legend.text=element_text(size=14), legend.position = 'top') +
  scale_fill_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(21, 21))

snpgdsClose(genofile)


### Xchr all variants (before any filtering)
vcf.fn = "chr8.gatk.called.vqsr.sv.placentas.vcf"
snpgdsVCF2GDS(vcf.fn, "chr8_all_variants.gds", method="biallelic.only")
snpgdsSummary("chr8_all_variants.gds")
genofile <- snpgdsOpen("chr8_all_variants.gds")
pca = snpgdsPCA(genofile, autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# Make a data frame
sample.id = pca$sample.id
chr8_df <- data.frame(sample.id = pca$sample.id,
                      pop = factor(pop_code)[match(pca$sample.id, sample.id)],
                      EV1 = pca$eigenvect[,1],    # the first eigenvector
                      EV2 = pca$eigenvect[,2],    # the second eigenvector
                      EV3 = pca$eigenvect[,3],
                      EV4 = pca$eigenvect[,4],
                      EV5 = pca$eigenvect[,5],
                      stringsAsFactors = FALSE)


p2 = ggplot(chr8_df, aes(x=EV1, y=EV2, fill=pop, shape=pop)) + 
  geom_point(size=3, colour="black") + 
  theme_bw()+
  labs(title="Chromosome 8", x="PC1", y="PC2") +
  theme(axis.text.x = element_text(size=14, colour="black"), axis.text.y = element_text(size=14, colour = "black"), axis.title.y = element_text(size=18, angle=0, vjust=0.5), axis.title.x = element_text(size=18), axis.line=element_line(colour = "black"),plot.title = element_text(hjust = 0.5, size=18), legend.title = element_blank(), panel.background = element_blank(), legend.text=element_text(size=14), legend.position = 'top') +
  scale_fill_manual(values=c("black", "red")) +
  scale_shape_manual(values=c(21, 21))

snpgdsClose(genofile)

png("/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_s3_pca.png", width=11, height=6, units = "in", res = 300)
ggarrange(p1, p2, ncol = 2, common.legend = T)
dev.off()

# pop = factor(pop_code)[match(pca$sample.id, sample.id)]
# 
# png("pca_plot.png", width=10, height=7, units = "in", res = 300)
# lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
# pairs(pca$eigenvect[,1:4], col=pop, labels=lbls)
# dev.off()
