rm(list = ls())
library(dplyr)
library(ggplot2)
## closest distance
YY1_processed <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/rev/processed_pgene/closest_with_L1/YY1_0.696_processed_pgene.txt')
Yy1_processed <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/rev/processed_pgene/closest_with_L1/Yy1_0.712_processed_pgene.txt')

YY1_unprocessed <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/rev/unprocessed_pgene/closest_with_L1/YY1_0.696_unprocessed_pgene.txt')
Yy1_unprocessed <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/rev/unprocessed_pgene/closest_with_L1/Yy1_0.712_unprocessed_pgene.txt')

YY1_lncRNA <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/rev/lncRNA/closest_with_L1/YY1_0.696_lncRNA.txt')
Yy1_lncRNA <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/rev/lncRNA/closest_with_L1/Yy1_0.712_lncRNA.txt')

YY1_pcg <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/rev/pcg/closest_with_L1/YY1_0.696_protein_coding.txt')
Yy1_pcg <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/rev/pcg/closest_with_L1/Yy1_0.712_protein_coding.txt')

df1 <- rbind(YY1_processed, Yy1_processed) %>% mutate(geneType = 'processed_pseudogene') %>%
  select(V22, geneType)
df2 <- rbind(YY1_unprocessed, Yy1_unprocessed) %>% mutate(geneType = 'unprocessed_pseudogene') %>%
  select(V22, geneType)
df3 <- rbind(YY1_lncRNA, Yy1_lncRNA) %>% mutate(geneType = 'lncRNA') %>%
  select(V22, geneType)
df4 <- rbind(YY1_pcg, Yy1_pcg) %>% mutate(geneType = 'protein_coding') %>%
  select(V22, geneType)

cols <- c('Protein-cod.' = 'black', 'lncRNA' = 'gray','Unprocessed' = '#edae49', 'Processed' = '#d1495b')

data_to_plot <- rbind(df1, df2, df3, df4)
data_to_plot$geneType <- factor(
  data_to_plot$geneType,
  levels = c("protein_coding", "lncRNA", "processed_pseudogene", "unprocessed_pseudogene"),
  labels = c("Protein-cod.", "lncRNA", "Processed", "Unprocessed"))

medians <- data_to_plot %>%
  group_by(geneType) %>%
  summarise(median_val = median(V22)/1000)

library(ggplot2)
library(ggdist)
library(gghalves)

my_comparisons <- list(
  c("Protein-cod.", "Processed"),
  c("lncRNA", "Processed"),
  c("Processed", "Unprocessed"))

p <- ggplot(data_to_plot, aes(x = geneType, y = V22/1000, fill = geneType)) +
  stat_halfeye(adjust = 0.4, width = 0.6, .width = 0, justification = -0.3, point_color = NA, slab_color = NA) +
  geom_boxplot(aes(color = geneType), width = 0.1, outlier.shape = NA, alpha = 0.6, position = position_nudge(x = 0)) +
  scale_color_manual(values = cols) + scale_fill_manual(values = cols) +
  labs(y = 'Distance to closest L1 (Kbp)', x = '', title = 'YY1-like motifs') +   
  scale_y_sqrt(breaks = pretty(c(0, max(data_to_plot$V22 / 1000)), n = 7)) + theme_bw() +
  ggpubr::stat_compare_means(comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif", , step.increase = 0.08, tip.length = 0.01) +
  theme(
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    axis.text = element_text(color = "black", size = 12),
    legend.position = 'none')
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/plots/closest_distance.pdf', p, width = 5, height = 4)

# overlap
rm(list = ls())
YY1_processed <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/rev/processed_pgene/overlap_with_L1/YY1_0.696_processed_pgene.txt')
Yy1_processed <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/rev/processed_pgene/overlap_with_L1/Yy1_0.712_processed_pgene.txt')
