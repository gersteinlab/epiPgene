rm(list = ls())
source('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/codes/PieDonutCustom.R')
library(dplyr)
library(tidyverse)
library(forcats)
cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')
setwd('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/intersection/HiC/processed/') # transcript or promoter region
files <- list.files(pattern = '*.txt')
df_list <- lapply(files, function(x) read.table(x, header = F, sep = '\t', stringsAsFactors = FALSE) %>% mutate(tissue = str_split(x, '.txt')[[1]][1]))
geneType <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T) %>% dplyr::select(gene_id, transcript_id, geneType) %>% distinct()

all <- bind_rows(df_list)
all <- all %>% dplyr::rename(gene_id = V10, transcript_id = V11)
all <- all %>% merge(geneType, by = c('gene_id', 'transcript_id'))

## PieDonut for pgenes
count <- all %>%
  select(gene_id, transcript_id, V1, V2, V3, V4, V5, V6, geneType) %>% distinct() %>% 
  mutate(interval_length = factor(paste((V6 - V5)/10^3, 'kbp'), levels = c('1 kbp', '2 kbp', '5 kbp', '10 kbp', '25 kbp'))) %>%
  group_by(geneType, interval_length) %>% summarise(n = n())

pdf('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/plots/pie_chart_pgene.pdf', width = 8, height = 8)
count %>% dplyr::filter(str_detect(geneType, 'pseudogene')) %>% 
  mutate(Pseudogene = str_remove(geneType, '_pseudogene') %>% str_to_title()) %>%
  PieDonutCustom(., aes(Pseudogene, interval_length, count = n), addDonutLabel = T, explodeDonut = T, explodePie = T, showPieName = F, addPieLabel = T, ratioByGroup = T, showRatioThreshold = 0.02,
           donutLabelSize = 5, pieLabelSize = 5, mainCol = c('#d1495b', '#00798c', '#edae49'))
dev.off()

library("viridis") 
p <- count %>% dplyr::filter(!str_detect(geneType, 'pseudogene')) %>%
  group_by(geneType) %>%
  mutate(prop = n / sum(n) * 100) %>%
  ggplot(aes(x = geneType, y = n, fill = interval_length)) +  
  geom_bar(stat = "identity", position = 'fill') +
  geom_text(aes(label = paste0(round(prop, 1), "%")), 
            position = position_fill(vjust = 0.5), 
            color = "white") +
  labs(x = '', y = 'Proportion', fill = 'Loop length') + theme_bw() +
  scale_x_discrete(labels = c('protein_coding' = 'Protein-coding', 'lncRNA' = 'lncRNA'), limits = c('protein_coding', 'lncRNA')) +
  scale_fill_viridis(discrete = TRUE, option = 'viridis', begin = 0, end = 0.6) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'right',
    axis.text = element_text(color="black", size = 12))

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/plots/barplot_PCG_lncRNA.pdf', p, width = 5, height = 5)


## LOOPS
library(ggpubr)
data_to_plot <- all %>%
  select(gene_id, transcript_id, V1, V2, V3, V4, V5, V6, V7, geneType) %>% distinct() %>%
  mutate(interval_length = factor(paste((V6 - V5)/10^3, 'kbp'), levels = c('1 kbp', '2 kbp', '5 kbp', '10 kbp', '25 kbp')))
data_to_plot[['geneType']] <-  factor(data_to_plot$geneType, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene'))

data_to_plot2 <- data_to_plot %>% # keep the pair with the max contact freq
  group_by(across(-V7)) %>% summarise(V7 = max(V7))

data_to_plot2 %>% ungroup() %>% rstatix::wilcox_test(V7~geneType, ref.group = 'processed_pseudogene', p.adjust.method = 'fdr')

p <- ggplot(data_to_plot2, aes(x = geneType, y = V7)) + 
  geom_violin(aes(fill = geneType), alpha = 1, trim = F)+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1),
               geom="pointrange", color="white", size = 0.1) +
#  facet_wrap(~interval_length, scales = 'free_x') +
  scale_x_discrete(labels = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed', 'Unitary')) + 
  scale_y_continuous(breaks = c(100, 200, 300, 400)) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  labs(x = '',  y = 'Observed contact frequency') +
  stat_compare_means(comparisons = list(c("processed_pseudogene", "protein_coding"), c("processed_pseudogene", "lncRNA"), c("processed_pseudogene", "unprocessed_pseudogene"), c("processed_pseudogene", "unitary_pseudogene")),
                     method = 'wilcox.test', label = 'p.signif', tip.length = 0.01, step.increase = 0.1,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'n.s.'))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'None',
    axis.text = element_text(color="black", size = 12))

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/plots/HiC_contact_freq_transcript.pdf', p, width = 6, height = 3)

library(rstatix)
library(dunn.test)
kw_contact <- data_to_plot2 %>% as.data.frame() %>% kruskal_test(V7 ~ geneType)
dunn.test(data_to_plot2$V7, data_to_plot2$geneType, method = 'bonferroni')
tbl <- data_to_plot2 %>% ungroup() %>% wilcox_test(V7~geneType, ref.group = 'processed_pseudogene', p.adjust.method = 'fdr')
write.table(tbl, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/tables/wilcoxon_dist.txt', sep = '\t', quote = F)

## draft
tmp <- all %>% dplyr::filter(geneType %in% c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene'))
write.table(tmp, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/tables/hic_coords_transcripts.txt', quote = F, sep = '\t', row.names = F, col.names = F)
