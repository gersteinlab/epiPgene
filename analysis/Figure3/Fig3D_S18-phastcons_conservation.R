rm(list = ls())
library(rstatix)
library(ggplot2)
library(tidyverse)
library(ggpubr)
conservation <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/phastCons/upstream_phastcons20way.bed')
colnames(conservation) <- c('chr', 'start', 'end', 'id', 'conservation')
mastertable <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/GCcontent/GRCh38_tss_upstream.v29.bed', header = F, sep = '\t')
colnames(mastertable) <- c('chr', 'start', 'end', 'gene_id', 'score', 'strand', 'transcript_id', 'geneType')
df <- mastertable %>% left_join(conservation, by = c('chr', 'start', 'end'))
cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')

df$geneType <- factor(df$geneType, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene'))

p <- df %>% filter(geneType != 'unitary_pseudogene') %>% ggplot(aes(x = geneType, y = conservation, fill = geneType, color = geneType)) + 
  geom_boxplot(alpha = 0.6, width = 0.4) + scale_y_sqrt() + labs(x = '', y = 'phastCons 20-way score', title = 'Upstream flanking region') + 
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  scale_x_discrete(labels = c(
    protein_coding = "Protein-cod.",
    lncRNA = "lncRNA",
    processed_pseudogene = "Processed",
    unprocessed_pseudogene = "Unprocessed")) +
  stat_compare_means(
    ref.group = 'processed_pseudogene',
    method = "wilcox.test") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'None',
    axis.text = element_text(color="black", size = 12))

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/plots/upstream_conservation_boxplot.pdf', p, width = 5.5, height = 4)


p <- df %>% filter(geneType != 'unitary_pseudogene') %>% ggplot(aes(x = conservation)) + 
  stat_ecdf(aes(color = factor(geneType, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene'))), linewidth = 0.6) + 
  scale_color_manual(values = cols,  label = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed', 'Unitary')) +
  labs(x = 'phastCons 20-way score', y = 'Cumulative frequency', color = 'Biotype', ) +
  theme_bw() + coord_fixed() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'right',
    axis.text = element_text(color="black", size = 12))
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/plots/upstream_conservation.pdf', p, width = 5, height = 5)

df$geneType <- factor(df$geneType, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene'))
p <- df %>% filter(geneType != 'unitary_pseudogene') %>% ggplot(aes(x = geneType, y = conservation, fill = geneType, color = geneType)) + 
  geom_boxplot(alpha = 0.6) +
  labs(x = '', y = 'Conservation',title = 'Upstream flanking region') +
  stat_compare_means(comparisons = list(c("processed_pseudogene", "protein_coding"), c("processed_pseudogene", "lncRNA"), c("processed_pseudogene", "unprocessed_pseudogene")),
                     method = 'wilcox.test', label = 'p.signif', tip.length = 0.01, step.increase = 0.05,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'n.s.'))) +
  scale_x_discrete(labels = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed'))+
  scale_y_continuous(breaks = seq(0,10,2) * 0.1) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'None',
    axis.text = element_text(color="black", size = 12))

tbl <- df %>% na.omit() %>% wilcox_test(conservation ~ geneType, p.adjust.method = 'fdr', ref.group = 'processed_pseudogene')
write.table(tbl, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/tables/upstream_conservation_wilcox_test.txt', sep = '\t', quote = F)
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/figS/upstream_conservation_boxplot.pdf', p, width = 5, height = 5)

## downstream
rm(list = ls())
conservation <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/phastCons/downstream_phastcons20way.bed')
colnames(conservation) <- c('chr', 'start', 'end', 'id', 'conservation')
mastertable <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/GCcontent/GRCh38_tss_downstream.v29.bed', header = F, sep = '\t')
colnames(mastertable) <- c('chr', 'start', 'end', 'gene_id', 'score', 'strand', 'transcript_id', 'geneType')
df <- mastertable %>% left_join(conservation, by = c('chr', 'start', 'end'))
cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')
p <- df %>% filter(geneType != 'unitary_pseudogene') %>% ggplot(aes(x = conservation)) + 
  stat_ecdf(aes(color = factor(geneType, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene'))), linewidth = 0.6) + 
  scale_color_manual(values = cols,  label = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed', 'Unitary')) +
  labs(x = 'phastCons 20-way score', y = 'Cumulative frequency', color = 'Biotype', title = 'Downstream flanking region') +
  theme_bw() + coord_fixed() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'right',
    axis.text = element_text(color="black", size = 12))

tbl <- df %>% na.omit() %>% wilcox_test(conservation ~ geneType, p.adjust.method = 'fdr', ref.group = 'processed_pseudogene')
write.table(tbl, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/tables/downstream_conservation_wilcox_test.txt', sep = '\t', quote = F)
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/plots/downstream_conservation.pdf', p, width = 5, height = 5)

df$geneType <- factor(df$geneType, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene'))
p <- df %>% filter(geneType != 'unitary_pseudogene') %>% ggplot(aes(x = geneType, y = conservation, fill = geneType, color = geneType)) + 
  geom_boxplot(alpha = 0.6) +
  labs(x = '', y = 'Conservation',title = 'Downstream flanking region') +
  stat_compare_means(comparisons = list(c("processed_pseudogene", "protein_coding"), c("processed_pseudogene", "lncRNA"), c("processed_pseudogene", "unprocessed_pseudogene")),
                     method = 'wilcox.test', label = 'p.signif', tip.length = 0.01, step.increase = 0.05,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'n.s.'))) +
  scale_x_discrete(labels = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed'))+
  scale_y_continuous(breaks = seq(0,10,2) * 0.1) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'None',
    axis.text = element_text(color="black", size = 12))
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/figS/downstream_conservation_boxplot.pdf', p, width = 5, height = 5)
