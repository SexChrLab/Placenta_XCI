library(tidyverse)
library(ggplot2)

setwd("/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/")
# -------------------------
# Compare gtex and placenta
# -------------------------
data = read.csv('categorize_genes_gtex_placenta_sorted.csv')
colnames(data) = c('Gene', 'Adult_Tissues', 'Placenta')

data_long <- gather(data, category, status, Placenta, Adult_Tissues,  factor_key=TRUE)

data_long$Gene <- as.character(data_long$Gene)
data_long$Gene <- factor(data_long$Gene, levels=unique(data_long$Gene))

png('/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_4A.png', width = 11, height = 4, units = "in", res = 1200)
ggplot(data_long, aes(category, Gene)) +
  geom_tile(aes(fill = status)) +
  labs(x="") +
  theme(axis.text.y = element_text(color = "black", size=16), legend.position = "none", axis.text.x=element_blank(), axis.title.x=element_text(size = 16), legend.text = element_text(size=18)) +
  scale_fill_manual(name = "", labels = c("Escape", "Inactivated", "Variable"), values = c("navyblue", "gold", "skyblue")) +
  coord_flip()
dev.off()

# legend only
p = ggplot(gtex_placenta_rmNA_sorted_long, aes(category, Gene)) +
  geom_tile(aes(fill = status)) +
  labs(x="") +
  theme(axis.text.y = element_text(color = "black", size=16), axis.text.x=element_blank(), axis.title.x=element_text(size = 16), legend.text = element_text(size=18)) +
  # scale_fill_manual(name = "", labels = c("Escape", "Inactivated", "Variable"), values = c("navyblue", "gold", "skyblue")) +
  scale_fill_manual(name = "", labels = c("Inactivated", "Escape", "Variable"), values = c("gold", "navyblue", "skyblue")) +
  coord_flip()
a = get_legend(p)
as_ggplot(a)


# Plot genes with opposite patterns between GTEx and Placenta
setwd('/scratch/tphung3/Placenta_XCI/gene_level/gtex_counts/genes_allele_balance/chrX/')

# Escape in placenta, inactivated in gtex: ['ACOT9', 'CASK', 'TSC22D3']
ACOT9_placenta = read.csv('ACOT9_allele_balance_for_placenta_samples.tsv', sep = '\t')
ACOT9_gtex = read.csv('ACOT9_allele_balance_for_gtex_samples.tsv', sep = '\t')

CASK_placenta = read.csv('CASK_allele_balance_for_placenta_samples.tsv', sep = '\t')
CASK_gtex = read.csv('CASK_allele_balance_for_gtex_samples.tsv', sep = '\t')

TSC22D3_placenta = read.csv('TSC22D3_allele_balance_for_placenta_samples.tsv', sep = '\t')
TSC22D3_gtex = read.csv('TSC22D3_allele_balance_for_gtex_samples.tsv', sep = '\t')

escape_placenta_inactivated_gtex = data.frame(
  allele_balance = c(ACOT9_placenta$allele_balance, ACOT9_gtex$allele_balance,
                     CASK_placenta$allele_balance, CASK_gtex$allele_balance,
                     TSC22D3_placenta$allele_balance, TSC22D3_gtex$allele_balance),
  labels = c(rep('Placenta', nrow(ACOT9_placenta)), rep('GTEx', nrow(ACOT9_gtex)),
             rep('Placenta', nrow(CASK_placenta)), rep('GTEx', nrow(CASK_gtex)),
             rep('Placenta', nrow(TSC22D3_placenta)), rep('GTEx', nrow(TSC22D3_gtex))),
  genes = c(rep('ACOT9', nrow(ACOT9_placenta) + nrow(ACOT9_gtex)),
            rep('CASK', nrow(CASK_placenta) + nrow(CASK_gtex)),
            rep('TSC22D3', nrow(TSC22D3_placenta) + nrow(TSC22D3_gtex))))
escape_placenta_inactivated_gtex$labels <- as.character(escape_placenta_inactivated_gtex$labels)
escape_placenta_inactivated_gtex$labels <- factor(escape_placenta_inactivated_gtex$labels, levels=unique(escape_placenta_inactivated_gtex$labels))

png('/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_4B.png', width = 4, height = 2, units = "in", res = 1200)
ggplot(escape_placenta_inactivated_gtex, aes(x=escape_placenta_inactivated_gtex$genes, y=escape_placenta_inactivated_gtex$allele_balance, color=labels)) +
  geom_boxplot(size=1.25) +
  theme_bw() +
  labs(x="", y="Allele balance", title = 'Placenta:Escape vs GTEx:Inactivated') +
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.title.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size=14), legend.position = "none") +
  coord_cartesian(ylim=c(0.5, 1)) +
  geom_hline(yintercept = 0.8, col="darkgray", linetype=2, size=1.5) +
  geom_hline(yintercept = 1, col="darkgray", linetype=2, size=1.5) +
  scale_color_manual(values=c("#C77CFF", "#00BFC4"))
dev.off()

# Obtain legends
p = ggplot(escape_placenta_inactivated_gtex, aes(x=escape_placenta_inactivated_gtex$genes, y=escape_placenta_inactivated_gtex$allele_balance, color=labels)) +
  geom_boxplot(size=1.25) +
  theme_bw() +
  labs(x="", y="Allele balance", title = 'Placenta:Escape vs GTEx:Inactivated') +
  theme(axis.text.x = element_text(size = 12), axis.title.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size=14), legend.text = element_text(size = 20), legend.title = element_blank()) +
  coord_cartesian(ylim=c(0.5, 1)) +
  geom_hline(yintercept = 0.8, col="darkgray", linetype=2, size=1.5) +
  geom_hline(yintercept = 1, col="darkgray", linetype=2, size=1.5) +
  scale_color_manual(values=c("#C77CFF", "#00BFC4"))
a = get_legend(p)
as_ggplot(a)

# Inactivated in placenta, escape in GTEx: ARSD, PRKX
ARSD_placenta = read.csv('ARSD_allele_balance_for_placenta_samples.tsv', sep = '\t')
ARSD_gtex = read.csv('ARSD_allele_balance_for_gtex_samples.tsv', sep = '\t')

PRKX_placenta = read.csv('PRKX_allele_balance_for_placenta_samples.tsv', sep = '\t')
PRKX_gtex = read.csv('PRKX_allele_balance_for_gtex_samples.tsv', sep = '\t')


inactivated_placenta_escape_gtex = data.frame(
  allele_balance = c(ARSD_placenta$allele_balance, ARSD_gtex$allele_balance,
                     PRKX_placenta$allele_balance, PRKX_gtex$allele_balance),
  labels = c(rep('Placenta', nrow(ARSD_placenta)), rep('GTEx', nrow(ARSD_gtex)),
            rep('Placenta', nrow(PRKX_placenta)), rep('GTEx', nrow(PRKX_gtex))),
  genes = c(rep('ARSD', nrow(ARSD_placenta) + nrow(ARSD_gtex)),
            rep('PRKX', nrow(PRKX_placenta) + nrow(PRKX_gtex))))
inactivated_placenta_escape_gtex$labels <- as.character(inactivated_placenta_escape_gtex$labels)
inactivated_placenta_escape_gtex$labels <- factor(inactivated_placenta_escape_gtex$labels, levels=unique(inactivated_placenta_escape_gtex$labels))

png('/scratch/tphung3/Placenta_XCI/manuscript_plots/plots/figure_4C.png', width = 4, height = 2, units = "in", res = 1200)
ggplot(inactivated_placenta_escape_gtex, aes(x=inactivated_placenta_escape_gtex$genes, y=inactivated_placenta_escape_gtex$allele_balance, color=labels)) +
  geom_boxplot(size=1.25) +
  theme_bw() +
  labs(x="", y="Allele balance", title = 'Placenta:Inactivated vs GTEx:Escape') +
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.title.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size=14), legend.position = "none") +
  coord_cartesian(ylim=c(0.5, 1)) +
  geom_hline(yintercept = 0.8, col="darkgray", linetype=2) +
  geom_hline(yintercept = 1, col="darkgray", linetype=2) +
  scale_color_manual(values=c("#C77CFF", "#00BFC4"))
dev.off()

# Variable escape in placenta, inactivated in GTEx: ['AL133500.1', 'ATG4A', 'ATRX', 'BCOR', 'CXorf56', 'EFNB1', 'GPC4', 'HCFC1', 'IL13RA1', 'IRAK1', 'LINC00632', 'MBTPS2', 'MSL3', 'PLS3', 'RBM41', 'RPS6KA3', 'SMS', 'TAB3', 'TSPAN6', 'USP11', 'VSIG4']
AL133500.1_placenta = read.csv('AL133500.1_allele_balance_for_placenta_samples.tsv', sep = '\t')
AL133500.1_gtex = read.csv('AL133500.1_allele_balance_for_gtex_samples.tsv', sep = '\t')
ATG4A_placenta = read.csv('ATG4A_allele_balance_for_placenta_samples.tsv', sep = '\t')
ATG4A_gtex = read.csv('ATG4A_allele_balance_for_gtex_samples.tsv', sep = '\t')
ATRX_placenta = read.csv('ATRX_allele_balance_for_placenta_samples.tsv', sep = '\t')
ATRX_gtex = read.csv('ATRX_allele_balance_for_gtex_samples.tsv', sep = '\t')
BCOR_placenta = read.csv('BCOR_allele_balance_for_placenta_samples.tsv', sep = '\t')
BCOR_gtex = read.csv('BCOR_allele_balance_for_gtex_samples.tsv', sep = '\t')
CXorf56_placenta = read.csv('CXorf56_allele_balance_for_placenta_samples.tsv', sep = '\t')
CXorf56_gtex = read.csv('CXorf56_allele_balance_for_gtex_samples.tsv', sep = '\t')
EFNB1_placenta = read.csv('EFNB1_allele_balance_for_placenta_samples.tsv', sep = '\t')
EFNB1_gtex = read.csv('EFNB1_allele_balance_for_gtex_samples.tsv', sep = '\t')
GPC4_placenta = read.csv('GPC4_allele_balance_for_placenta_samples.tsv', sep = '\t')
GPC4_gtex = read.csv('GPC4_allele_balance_for_gtex_samples.tsv', sep = '\t')
HCFC1_placenta = read.csv('HCFC1_allele_balance_for_placenta_samples.tsv', sep = '\t')
HCFC1_gtex = read.csv('HCFC1_allele_balance_for_gtex_samples.tsv', sep = '\t')
IL13RA1_placenta = read.csv('IL13RA1_allele_balance_for_placenta_samples.tsv', sep = '\t')
IL13RA1_gtex = read.csv('IL13RA1_allele_balance_for_gtex_samples.tsv', sep = '\t')
IRAK1_placenta = read.csv('IRAK1_allele_balance_for_placenta_samples.tsv', sep = '\t')
IRAK1_gtex = read.csv('IRAK1_allele_balance_for_gtex_samples.tsv', sep = '\t')
LINC00632_placenta = read.csv('LINC00632_allele_balance_for_placenta_samples.tsv', sep = '\t')
LINC00632_gtex = read.csv('LINC00632_allele_balance_for_gtex_samples.tsv', sep = '\t')
MBTPS2_placenta = read.csv('MBTPS2_allele_balance_for_placenta_samples.tsv', sep = '\t')
MBTPS2_gtex = read.csv('MBTPS2_allele_balance_for_gtex_samples.tsv', sep = '\t')
MSL3_placenta = read.csv('MSL3_allele_balance_for_placenta_samples.tsv', sep = '\t')
MSL3_gtex = read.csv('MSL3_allele_balance_for_gtex_samples.tsv', sep = '\t')
PLS3_placenta = read.csv('PLS3_allele_balance_for_placenta_samples.tsv', sep = '\t')
PLS3_gtex = read.csv('PLS3_allele_balance_for_gtex_samples.tsv', sep = '\t')
RBM41_placenta = read.csv('RBM41_allele_balance_for_placenta_samples.tsv', sep = '\t')
RBM41_gtex = read.csv('RBM41_allele_balance_for_gtex_samples.tsv', sep = '\t')
RPS6KA3_placenta = read.csv('RPS6KA3_allele_balance_for_placenta_samples.tsv', sep = '\t')
RPS6KA3_gtex = read.csv('RPS6KA3_allele_balance_for_gtex_samples.tsv', sep = '\t')
SMS_placenta = read.csv('SMS_allele_balance_for_placenta_samples.tsv', sep = '\t')
SMS_gtex = read.csv('SMS_allele_balance_for_gtex_samples.tsv', sep = '\t')
TAB3_placenta = read.csv('TAB3_allele_balance_for_placenta_samples.tsv', sep = '\t')
TAB3_gtex = read.csv('TAB3_allele_balance_for_gtex_samples.tsv', sep = '\t')
TSPAN6_placenta = read.csv('TSPAN6_allele_balance_for_placenta_samples.tsv', sep = '\t')
TSPAN6_gtex = read.csv('TSPAN6_allele_balance_for_gtex_samples.tsv', sep = '\t')
USP11_placenta = read.csv('USP11_allele_balance_for_placenta_samples.tsv', sep = '\t')
USP11_gtex = read.csv('USP11_allele_balance_for_gtex_samples.tsv', sep = '\t')
VSIG4_placenta = read.csv('VSIG4_allele_balance_for_placenta_samples.tsv', sep = '\t')
VSIG4_gtex = read.csv('VSIG4_allele_balance_for_gtex_samples.tsv', sep = '\t')

variable_placenta_inactivated_gtex = data.frame(
  allele_balance = c(AL133500.1_placenta$allele_balance, AL133500.1_gtex$allele_balance,
                     ATG4A_placenta$allele_balance, ATG4A_gtex$allele_balance,
                     ATRX_placenta$allele_balance, ATRX_gtex$allele_balance,
                     BCOR_placenta$allele_balance, BCOR_gtex$allele_balance,
                     CXorf56_placenta$allele_balance, CXorf56_gtex$allele_balance,
                     EFNB1_placenta$allele_balance, EFNB1_gtex$allele_balance,
                     GPC4_placenta$allele_balance, GPC4_gtex$allele_balance,
                     HCFC1_placenta$allele_balance, HCFC1_gtex$allele_balance,
                     IL13RA1_placenta$allele_balance, IL13RA1_gtex$allele_balance,
                     IRAK1_placenta$allele_balance, IRAK1_gtex$allele_balance,
                     LINC00632_placenta$allele_balance, LINC00632_gtex$allele_balance,
                     MBTPS2_placenta$allele_balance, MBTPS2_gtex$allele_balance,
                     MSL3_placenta$allele_balance, MSL3_gtex$allele_balance,
                     PLS3_placenta$allele_balance, PLS3_gtex$allele_balance,
                     RBM41_placenta$allele_balance, RBM41_gtex$allele_balance,
                     RPS6KA3_placenta$allele_balance, RPS6KA3_gtex$allele_balance,
                     SMS_placenta$allele_balance, SMS_gtex$allele_balance,
                     TAB3_placenta$allele_balance, TAB3_gtex$allele_balance,
                     TSPAN6_placenta$allele_balance, TSPAN6_gtex$allele_balance,
                     USP11_placenta$allele_balance, USP11_gtex$allele_balance,
                     VSIG4_placenta$allele_balance, VSIG4_gtex$allele_balance),
  labels = c(rep('Placenta', nrow(AL133500.1_placenta)), rep('GTEx', nrow(AL133500.1_gtex)),
             rep('Placenta', nrow(ATG4A_placenta)), rep('GTEx', nrow(ATG4A_gtex)),
             rep('Placenta', nrow(ATRX_placenta)), rep('GTEx', nrow(ATRX_gtex)),
             rep('Placenta', nrow(BCOR_placenta)), rep('GTEx', nrow(BCOR_gtex)),
             rep('Placenta', nrow(CXorf56_placenta)), rep('GTEx', nrow(CXorf56_gtex)),
             rep('Placenta', nrow(EFNB1_placenta)), rep('GTEx', nrow(EFNB1_gtex)),
             rep('Placenta', nrow(GPC4_placenta)), rep('GTEx', nrow(GPC4_gtex)),
             rep('Placenta', nrow(HCFC1_placenta)), rep('GTEx', nrow(HCFC1_gtex)),
             rep('Placenta', nrow(IL13RA1_placenta)), rep('GTEx', nrow(IL13RA1_gtex)),
             rep('Placenta', nrow(IRAK1_placenta)), rep('GTEx', nrow(IRAK1_gtex)),
             rep('Placenta', nrow(LINC00632_placenta)), rep('GTEx', nrow(LINC00632_gtex)),
             rep('Placenta', nrow(MBTPS2_placenta)), rep('GTEx', nrow(MBTPS2_gtex)),
             rep('Placenta', nrow(MSL3_placenta)), rep('GTEx', nrow(MSL3_gtex)),
             rep('Placenta', nrow(PLS3_placenta)), rep('GTEx', nrow(PLS3_gtex)),
             rep('Placenta', nrow(RBM41_placenta)), rep('GTEx', nrow(RBM41_gtex)),
             rep('Placenta', nrow(RPS6KA3_placenta)), rep('GTEx', nrow(RPS6KA3_gtex)),
             rep('Placenta', nrow(SMS_placenta)), rep('GTEx', nrow(SMS_gtex)),
             rep('Placenta', nrow(TAB3_placenta)), rep('GTEx', nrow(TAB3_gtex)),
             rep('Placenta', nrow(TSPAN6_placenta)), rep('GTEx', nrow(TSPAN6_gtex)),
             rep('Placenta', nrow(USP11_placenta)), rep('GTEx', nrow(USP11_gtex)),
             rep('Placenta', nrow(VSIG4_placenta)), rep('GTEx', nrow(VSIG4_gtex))),
  genes = c(rep('AL133500.1', nrow(AL133500.1_placenta) + nrow(AL133500.1_gtex)),
            rep('ATG4A', nrow(ATG4A_placenta) + nrow(ATG4A_gtex)),
            rep('ATRX', nrow(ATRX_placenta) + nrow(ATRX_gtex)),
            rep('BCOR', nrow(BCOR_placenta) + nrow(BCOR_gtex)),
            rep('CXorf56', nrow(CXorf56_placenta) + nrow(CXorf56_gtex)),
            rep('EFNB1', nrow(EFNB1_placenta) + nrow(EFNB1_gtex)),
            rep('GPC4', nrow(GPC4_placenta) + nrow(GPC4_gtex)),
            rep('HCFC1', nrow(HCFC1_placenta) + nrow(HCFC1_gtex)),
            rep('IL13RA1', nrow(IL13RA1_placenta) + nrow(IL13RA1_gtex)),
            rep('IRAK1', nrow(IRAK1_placenta) + nrow(IRAK1_gtex)),
            rep('LINC00632', nrow(LINC00632_placenta) + nrow(LINC00632_gtex)),
            rep('MBTPS2', nrow(MBTPS2_placenta) + nrow(MBTPS2_gtex)),
            rep('MSL3', nrow(MSL3_placenta) + nrow(MSL3_gtex)),
            rep('PLS3', nrow(PLS3_placenta) + nrow(PLS3_gtex)),
            rep('RBM41', nrow(RBM41_placenta) + nrow(RBM41_gtex)),
            rep('RPS6KA3', nrow(RPS6KA3_placenta) + nrow(RPS6KA3_gtex)),
            rep('SMS', nrow(SMS_placenta) + nrow(SMS_gtex)),
            rep('TAB3', nrow(TAB3_placenta) + nrow(TAB3_gtex)),
            rep('TSPAN6', nrow(TSPAN6_placenta) + nrow(TSPAN6_gtex)),
            rep('USP11', nrow(USP11_placenta) + nrow(USP11_gtex)),
            rep('VSIG4', nrow(VSIG4_placenta) + nrow(VSIG4_gtex))))

variable_placenta_inactivated_gtex$labels <- as.character(variable_placenta_inactivated_gtex$labels)
variable_placenta_inactivated_gtex$labels <- factor(variable_placenta_inactivated_gtex$labels, levels=unique(variable_placenta_inactivated_gtex$labels))


ggplot(variable_placenta_inactivated_gtex, aes(x=variable_placenta_inactivated_gtex$genes, y=variable_placenta_inactivated_gtex$allele_balance, color=labels)) +
  geom_boxplot(size=1.75)+
  theme_bw() +
  labs(x="", y="Allele balance", title = 'Placenta:Variable vs GTEx:Inactivated') +
  theme(axis.text.x = element_text(size = 30, colour = "black", angle=30), axis.title.y = element_text(size = 36), plot.title = element_text(hjust = 0.5, size=46), legend.position = "none") +
  coord_cartesian(ylim=c(0.5, 1)) +
  geom_hline(yintercept = 0.8, col="darkgray", linetype=2) +
  geom_hline(yintercept = 1, col="darkgray", linetype=2) +
  scale_color_manual(values=c("#C77CFF", "#00BFC4"))


# Variable escape in placenta, escape in GTEx: ['CD99', 'PIN4', 'TXLNG', 'ZFX', 'ZRSR2']
CD99_placenta = read.csv('CD99_allele_balance_for_placenta_samples.tsv', sep = '\t')
CD99_gtex = read.csv('CD99_allele_balance_for_gtex_samples.tsv', sep = '\t')
PIN4_placenta = read.csv('PIN4_allele_balance_for_placenta_samples.tsv', sep = '\t')
PIN4_gtex = read.csv('PIN4_allele_balance_for_gtex_samples.tsv', sep = '\t')
TXLNG_placenta = read.csv('TXLNG_allele_balance_for_placenta_samples.tsv', sep = '\t')
TXLNG_gtex = read.csv('TXLNG_allele_balance_for_gtex_samples.tsv', sep = '\t')
ZFX_placenta = read.csv('ZFX_allele_balance_for_placenta_samples.tsv', sep = '\t')
ZFX_gtex = read.csv('ZFX_allele_balance_for_gtex_samples.tsv', sep = '\t')
ZRSR2_placenta = read.csv('ZRSR2_allele_balance_for_placenta_samples.tsv', sep = '\t')
ZRSR2_gtex = read.csv('ZRSR2_allele_balance_for_gtex_samples.tsv', sep = '\t')

variable_placenta_escape_gtex = data.frame(
  allele_balance = c(CD99_placenta$allele_balance, CD99_gtex$allele_balance,
                     PIN4_placenta$allele_balance, PIN4_gtex$allele_balance,
                     TXLNG_placenta$allele_balance, TXLNG_gtex$allele_balance,
                     ZFX_placenta$allele_balance, ZFX_gtex$allele_balance,
                     ZRSR2_placenta$allele_balance, ZRSR2_gtex$allele_balance),
  labels = c(rep('Placenta', nrow(CD99_placenta)), rep('GTEx', nrow(CD99_gtex)),
             rep('Placenta', nrow(PIN4_placenta)), rep('GTEx', nrow(PIN4_gtex)),
             rep('Placenta', nrow(TXLNG_placenta)), rep('GTEx', nrow(TXLNG_gtex)),
             rep('Placenta', nrow(ZFX_placenta)), rep('GTEx', nrow(ZFX_gtex)),
             rep('Placenta', nrow(ZRSR2_placenta)), rep('GTEx', nrow(ZRSR2_gtex))),
  genes = c(rep('CD99', nrow(CD99_placenta) + nrow(CD99_gtex)),
            rep('PIN4', nrow(PIN4_placenta) + nrow(PIN4_gtex)),
            rep('TXLNG', nrow(TXLNG_placenta) + nrow(TXLNG_gtex)),
            rep('ZFX', nrow(ZFX_placenta) + nrow(ZFX_gtex)),
            rep('ZRSR2', nrow(ZRSR2_placenta) + nrow(ZRSR2_gtex))))

variable_placenta_escape_gtex$labels <- as.character(variable_placenta_escape_gtex$labels)
variable_placenta_escape_gtex$labels <- factor(variable_placenta_escape_gtex$labels, levels=unique(variable_placenta_escape_gtex$labels))


ggplot(variable_placenta_escape_gtex, aes(x=variable_placenta_escape_gtex$genes, y=variable_placenta_escape_gtex$allele_balance, color=labels)) +
  geom_boxplot(size=1.25) +
  theme_bw() +
  labs(x="", y="Allele balance", title = 'Placenta:Variable vs GTEx:Escape') +
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.title.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size=14), legend.position = "none") +
  coord_cartesian(ylim=c(0.5, 1)) +
  geom_hline(yintercept = 0.8, col="darkgray", linetype=2) +
  geom_hline(yintercept = 1, col="darkgray", linetype=2) +
  scale_color_manual(values=c("#C77CFF", "#00BFC4"))


# Inactivated in placenta, variable escape in GTEx: ['ALG13', 'APLN', 'BGN', 'CAPN6', 'CHM', 'FHL1', 'GPC3', 'GPRASP2', 'HAUS7', 'NLGN4X', 'OFD1', 'PIR', 'PRRG1', 'SEPT6', 'SH3KBP1', 'STARD8', 'USP9X']
ALG13_placenta = read.csv('ALG13_allele_balance_for_placenta_samples.tsv', sep = '\t')
ALG13_gtex = read.csv('ALG13_allele_balance_for_gtex_samples.tsv', sep = '\t')
APLN_placenta = read.csv('APLN_allele_balance_for_placenta_samples.tsv', sep = '\t')
APLN_gtex = read.csv('APLN_allele_balance_for_gtex_samples.tsv', sep = '\t')
BGN_placenta = read.csv('BGN_allele_balance_for_placenta_samples.tsv', sep = '\t')
BGN_gtex = read.csv('BGN_allele_balance_for_gtex_samples.tsv', sep = '\t')
CAPN6_placenta = read.csv('CAPN6_allele_balance_for_placenta_samples.tsv', sep = '\t')
CAPN6_gtex = read.csv('CAPN6_allele_balance_for_gtex_samples.tsv', sep = '\t')
CHM_placenta = read.csv('CHM_allele_balance_for_placenta_samples.tsv', sep = '\t')
CHM_gtex = read.csv('CHM_allele_balance_for_gtex_samples.tsv', sep = '\t')
FHL1_placenta = read.csv('FHL1_allele_balance_for_placenta_samples.tsv', sep = '\t')
FHL1_gtex = read.csv('FHL1_allele_balance_for_gtex_samples.tsv', sep = '\t')
GPC3_placenta = read.csv('GPC3_allele_balance_for_placenta_samples.tsv', sep = '\t')
GPC3_gtex = read.csv('GPC3_allele_balance_for_gtex_samples.tsv', sep = '\t')
GPRASP2_placenta = read.csv('GPRASP2_allele_balance_for_placenta_samples.tsv', sep = '\t')
GPRASP2_gtex = read.csv('GPRASP2_allele_balance_for_gtex_samples.tsv', sep = '\t')
HAUS7_placenta = read.csv('HAUS7_allele_balance_for_placenta_samples.tsv', sep = '\t')
HAUS7_gtex = read.csv('HAUS7_allele_balance_for_gtex_samples.tsv', sep = '\t')
NLGN4X_placenta = read.csv('NLGN4X_allele_balance_for_placenta_samples.tsv', sep = '\t')
NLGN4X_gtex = read.csv('NLGN4X_allele_balance_for_gtex_samples.tsv', sep = '\t')
OFD1_placenta = read.csv('OFD1_allele_balance_for_placenta_samples.tsv', sep = '\t')
OFD1_gtex = read.csv('OFD1_allele_balance_for_gtex_samples.tsv', sep = '\t')
PIR_placenta = read.csv('PIR_allele_balance_for_placenta_samples.tsv', sep = '\t')
PIR_gtex = read.csv('PIR_allele_balance_for_gtex_samples.tsv', sep = '\t')
PRRG1_placenta = read.csv('PRRG1_allele_balance_for_placenta_samples.tsv', sep = '\t')
PRRG1_gtex = read.csv('PRRG1_allele_balance_for_gtex_samples.tsv', sep = '\t')
SEPT6_placenta = read.csv('SEPT6_allele_balance_for_placenta_samples.tsv', sep = '\t')
SEPT6_gtex = read.csv('SEPT6_allele_balance_for_gtex_samples.tsv', sep = '\t')
SH3KBP1_placenta = read.csv('SH3KBP1_allele_balance_for_placenta_samples.tsv', sep = '\t')
SH3KBP1_gtex = read.csv('SH3KBP1_allele_balance_for_gtex_samples.tsv', sep = '\t')
STARD8_placenta = read.csv('STARD8_allele_balance_for_placenta_samples.tsv', sep = '\t')
STARD8_gtex = read.csv('STARD8_allele_balance_for_gtex_samples.tsv', sep = '\t')
USP9X_placenta = read.csv('USP9X_allele_balance_for_placenta_samples.tsv', sep = '\t')
USP9X_gtex = read.csv('USP9X_allele_balance_for_gtex_samples.tsv', sep = '\t')

inactivated_placenta_variable_gtex = data.frame(
  allele_balance = c(ALG13_placenta$allele_balance, ALG13_gtex$allele_balance,
                     APLN_placenta$allele_balance, APLN_gtex$allele_balance,
                     BGN_placenta$allele_balance, BGN_gtex$allele_balance,
                     CAPN6_placenta$allele_balance, CAPN6_gtex$allele_balance,
                     CHM_placenta$allele_balance, CHM_gtex$allele_balance,
                     FHL1_placenta$allele_balance, FHL1_gtex$allele_balance,
                     GPC3_placenta$allele_balance, GPC3_gtex$allele_balance,
                     GPRASP2_placenta$allele_balance, GPRASP2_gtex$allele_balance,
                     HAUS7_placenta$allele_balance, HAUS7_gtex$allele_balance,
                     NLGN4X_placenta$allele_balance, NLGN4X_gtex$allele_balance,
                     OFD1_placenta$allele_balance, OFD1_gtex$allele_balance,
                     PIR_placenta$allele_balance, PIR_gtex$allele_balance,
                     PRRG1_placenta$allele_balance, PRRG1_gtex$allele_balance,
                     SEPT6_placenta$allele_balance, SEPT6_gtex$allele_balance,
                     SH3KBP1_placenta$allele_balance, SH3KBP1_gtex$allele_balance,
                     STARD8_placenta$allele_balance, STARD8_gtex$allele_balance,
                     USP9X_placenta$allele_balance, USP9X_gtex$allele_balance),
  labels = c(rep('Placenta', nrow(ALG13_placenta)), rep('GTEx', nrow(ALG13_gtex)),
             rep('Placenta', nrow(APLN_placenta)), rep('GTEx', nrow(APLN_gtex)),
             rep('Placenta', nrow(BGN_placenta)), rep('GTEx', nrow(BGN_gtex)),
             rep('Placenta', nrow(CAPN6_placenta)), rep('GTEx', nrow(CAPN6_gtex)),
             rep('Placenta', nrow(CHM_placenta)), rep('GTEx', nrow(CHM_gtex)),
             rep('Placenta', nrow(FHL1_placenta)), rep('GTEx', nrow(FHL1_gtex)),
             rep('Placenta', nrow(GPC3_placenta)), rep('GTEx', nrow(GPC3_gtex)),
             rep('Placenta', nrow(GPRASP2_placenta)), rep('GTEx', nrow(GPRASP2_gtex)),
             rep('Placenta', nrow(HAUS7_placenta)), rep('GTEx', nrow(HAUS7_gtex)),
             rep('Placenta', nrow(NLGN4X_placenta)), rep('GTEx', nrow(NLGN4X_gtex)),
             rep('Placenta', nrow(OFD1_placenta)), rep('GTEx', nrow(OFD1_gtex)),
             rep('Placenta', nrow(PIR_placenta)), rep('GTEx', nrow(PIR_gtex)),
             rep('Placenta', nrow(PRRG1_placenta)), rep('GTEx', nrow(PRRG1_gtex)),
             rep('Placenta', nrow(SEPT6_placenta)), rep('GTEx', nrow(SEPT6_gtex)),
             rep('Placenta', nrow(SH3KBP1_placenta)), rep('GTEx', nrow(SH3KBP1_gtex)),
             rep('Placenta', nrow(STARD8_placenta)), rep('GTEx', nrow(STARD8_gtex)),
             rep('Placenta', nrow(USP9X_placenta)), rep('GTEx', nrow(USP9X_gtex))),
  genes = c(rep('ALG13', nrow(ALG13_placenta) + nrow(ALG13_gtex)),
            rep('APLN', nrow(APLN_placenta) + nrow(APLN_gtex)),
            rep('BGN', nrow(BGN_placenta) + nrow(BGN_gtex)),
            rep('CAPN6', nrow(CAPN6_placenta) + nrow(CAPN6_gtex)),
            rep('CHM', nrow(CHM_placenta) + nrow(CHM_gtex)),
            rep('FHL1', nrow(FHL1_placenta) + nrow(FHL1_gtex)),
            rep('GPC3', nrow(GPC3_placenta) + nrow(GPC3_gtex)),
            rep('GPRASP2', nrow(GPRASP2_placenta) + nrow(GPRASP2_gtex)),
            rep('HAUS7', nrow(HAUS7_placenta) + nrow(HAUS7_gtex)),
            rep('NLGN4X', nrow(NLGN4X_placenta) + nrow(NLGN4X_gtex)),
            rep('OFD1', nrow(OFD1_placenta) + nrow(OFD1_gtex)),
            rep('PIR', nrow(PIR_placenta) + nrow(PIR_gtex)),
            rep('PRRG1', nrow(PRRG1_placenta) + nrow(PRRG1_gtex)),
            rep('SEPT6', nrow(SEPT6_placenta) + nrow(SEPT6_gtex)),
            rep('SH3KBP1', nrow(SH3KBP1_placenta) + nrow(SH3KBP1_gtex)),
            rep('STARD8', nrow(STARD8_placenta) + nrow(STARD8_gtex)),
            rep('USP9X', nrow(USP9X_placenta) + nrow(USP9X_gtex))))

inactivated_placenta_variable_gtex$labels <- as.character(inactivated_placenta_variable_gtex$labels)
inactivated_placenta_variable_gtex$labels <- factor(inactivated_placenta_variable_gtex$labels, levels=unique(inactivated_placenta_variable_gtex$labels))


ggplot(inactivated_placenta_variable_gtex, aes(x=inactivated_placenta_variable_gtex$genes, y=inactivated_placenta_variable_gtex$allele_balance, color=labels)) +
  geom_boxplot(size=1.25) +
  theme_bw() +
  labs(x="", y="Allele balance", title = 'Placenta:Inactivated vs GTEx:Variable') +
  theme(axis.text.x = element_text(size = 30, colour = "black", angle=30), axis.title.y = element_text(size = 36), plot.title = element_text(hjust = 0.5, size=46), legend.position = "none") +
  coord_cartesian(ylim=c(0.5, 1)) +
  geom_hline(yintercept = 0.8, col="darkgray", linetype=2) +
  geom_hline(yintercept = 1, col="darkgray", linetype=2) +
  scale_color_manual(values=c("#C77CFF", "#00BFC4"))


# Escape in placenta, variable escape in GTEx: ['ARHGAP4', 'CXorf36', 'FLNA', 'SMC1A', 'UBA1']
ARHGAP4_placenta = read.csv('ARHGAP4_allele_balance_for_placenta_samples.tsv', sep = '\t')
ARHGAP4_gtex = read.csv('ARHGAP4_allele_balance_for_gtex_samples.tsv', sep = '\t')
CXorf36_placenta = read.csv('CXorf36_allele_balance_for_placenta_samples.tsv', sep = '\t')
CXorf36_gtex = read.csv('CXorf36_allele_balance_for_gtex_samples.tsv', sep = '\t')
FLNA_placenta = read.csv('FLNA_allele_balance_for_placenta_samples.tsv', sep = '\t')
FLNA_gtex = read.csv('FLNA_allele_balance_for_gtex_samples.tsv', sep = '\t')
SMC1A_placenta = read.csv('SMC1A_allele_balance_for_placenta_samples.tsv', sep = '\t')
SMC1A_gtex = read.csv('SMC1A_allele_balance_for_gtex_samples.tsv', sep = '\t')
UBA1_placenta = read.csv('UBA1_allele_balance_for_placenta_samples.tsv', sep = '\t')
UBA1_gtex = read.csv('UBA1_allele_balance_for_gtex_samples.tsv', sep = '\t')

escape_placenta_variable_gtex = data.frame(
  allele_balance = c(ARHGAP4_placenta$allele_balance, ARHGAP4_gtex$allele_balance,
                     CXorf36_placenta$allele_balance, CXorf36_gtex$allele_balance,
                     FLNA_placenta$allele_balance, FLNA_gtex$allele_balance,
                     SMC1A_placenta$allele_balance, SMC1A_gtex$allele_balance,
                     UBA1_placenta$allele_balance, UBA1_gtex$allele_balance),
  labels = c(rep('Placenta', nrow(ARHGAP4_placenta)), rep('GTEx', nrow(ARHGAP4_gtex)),
             rep('Placenta', nrow(CXorf36_placenta)), rep('GTEx', nrow(CXorf36_gtex)),
             rep('Placenta', nrow(FLNA_placenta)), rep('GTEx', nrow(FLNA_gtex)),
             rep('Placenta', nrow(SMC1A_placenta)), rep('GTEx', nrow(SMC1A_gtex)),
             rep('Placenta', nrow(UBA1_placenta)), rep('GTEx', nrow(UBA1_gtex))),
  genes = c(rep('ARHGAP4', nrow(ARHGAP4_placenta) + nrow(ARHGAP4_gtex)),
            rep('CXorf36', nrow(CXorf36_placenta) + nrow(CXorf36_gtex)),
            rep('FLNA', nrow(FLNA_placenta) + nrow(FLNA_gtex)),
            rep('SMC1A', nrow(SMC1A_placenta) + nrow(SMC1A_gtex)),
            rep('UBA1', nrow(UBA1_placenta) + nrow(UBA1_gtex))))

escape_placenta_variable_gtex$labels <- as.character(escape_placenta_variable_gtex$labels)
escape_placenta_variable_gtex$labels <- factor(escape_placenta_variable_gtex$labels, levels=unique(escape_placenta_variable_gtex$labels))


ggplot(escape_placenta_variable_gtex, aes(x=escape_placenta_variable_gtex$genes, y=escape_placenta_variable_gtex$allele_balance, color=labels)) +
  geom_boxplot(size=1.25) +
  theme_bw() +
  labs(x="", y="Allele balance", title = 'Placenta:Escape vs GTEx:Variable') +
  theme(axis.text.x = element_text(size = 12, colour = "black"), axis.title.y = element_text(size = 12), plot.title = element_text(hjust = 0.5, size=14), legend.position = "none") +
  coord_cartesian(ylim=c(0.5, 1)) +
  geom_hline(yintercept = 0.8, col="darkgray", linetype=2) +
  geom_hline(yintercept = 1, col="darkgray", linetype=2) +
  scale_color_manual(values=c("#C77CFF", "#00BFC4"))


# -------------------------------------------
# Compare between placenta and prior research
# -------------------------------------------
# sorted:
placenta_prior_sorted = read.csv('c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/E_escape_genes/gtex_counts/fulltermplacenta_vs_priorresearch_sorted.csv', sep = ',')
placenta_prior_sorted_long <- gather(placenta_prior_sorted, category, status, Placenta, Prior_research,  factor_key=TRUE)

placenta_prior_sorted_long$Gene <- as.character(placenta_prior_sorted_long$Gene)
placenta_prior_sorted_long$Gene <- factor(placenta_prior_sorted_long$Gene, levels=unique(placenta_prior_sorted_long$Gene))

png('c://Users/tuyen/Documents/postdoc_asu/projects/PlacentaSexDiff/E_escape_genes/gtex_counts/plots/fulltermplacenta_vs_priorresearch.png', width = 6, height = 11, units = "in", res = 600)
ggplot(placenta_prior_sorted_long, aes(category, Gene)) +
  geom_tile(aes(fill = status)) +
  labs(x="", title="Placenta vs Prior studies") +
  theme(axis.text.y = element_blank(), legend.position = "top", axis.text.x=element_text(color = "black", size=16), axis.title.y=element_text(size = 16), legend.text = element_text(size=14), plot.title = element_text(hjust = 0.5, size=18)) +
  scale_fill_manual(name = "", labels = c("Escape", "Inactivated", "Variable"), values = c("navyblue", "gold", "skyblue"))
dev.off()

