rm(list = ls())
library(dplyr)
library(tidyverse)
library(ggpubr)
all <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
data_to_plot <- all %>% 
  filter(RAMPAGE == "True") %>%
  dplyr::select(geneType, transcript_id, tissue, avgTPM) %>% distinct() %>% mutate(expr = log1p(avgTPM))

cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b')

# stats
library(rstatix)
stat.test <- data_to_plot %>%
  filter(geneType != "unitary_pseudogene") %>%
  group_by(tissue) %>%
  wilcox_test(expr ~ geneType, ref.group = "processed_pseudogene") %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>% 
  mutate(y.position = 10)

data_to_plot$geneType <- factor(data_to_plot$geneType, levels = c('protein_coding','lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene'))

p <- data_to_plot %>% filter(geneType != 'unitary_pseudogene') %>%
  ggplot(aes(x = geneType, y = expr)) + geom_violin(aes(fill=geneType), size = 0.3, alpha = 0.6, trim = F, scale = "width") +
  geom_boxplot(width = 0.05, outlier.shape = NA, size = 0.3) + facet_wrap(~tissue) +
  scale_y_continuous(breaks = seq(1, 12, 2)) +
  coord_cartesian(ylim = c(-1, 12)) +
  scale_fill_manual(values = cols, labels = c(
    "lncRNA" = "lncRNA",
    "processed_pseudogene" = "Processed",
    "protein_coding" = "Protein-coding",
    "unprocessed_pseudogene" = "Unprocessed")) +
  labs(x = '', y = 'Log-transformed transcript expression', fill = 'Gene type') + theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color="black", size = 12),
        legend.position = "bottom")

p <- p + stat_pvalue_manual(
  stat.test %>% filter(p.adj < 0.05),
  label = "p.adj.signif",
  tip.length = 0.05,
  step.increase = 0)

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/tissueByTissue/plots/expression.pdf', p, width = 11, height = 8.5, units = 'in')
