df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
cols_1 <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'pseudogene' = '#BEB6DD')
cols_2 <- c('protein_coding' = 'black', 'lncRNA' = 'gray','unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')
library(dplyr)
library(stringr)
library(reshape2)
library(tidyverse)
long_data <- df %>% filter(RAMPAGE == "True") %>%
  mutate(geneType = if_else(str_detect(geneType, 'pseudogene'), 'pseudogene', geneType)) %>%
  mutate(tss_length = end - start) %>%
  pivot_longer(cols = c(H3K27ac, H3K27me3, H3K4me1, H3K4me3, H3K9me3, H3K36me3, DNase, ATAC), names_to = 'assay', values_to = 'coverage')

results <- long_data %>% na.omit() %>%
  group_by(geneType, assay, tissue) %>%
  summarize(
    avg_prop = mean(coverage / tss_length),
    std_prop = sd(coverage / tss_length)) %>%
  select(geneType, assay, tissue, avg_prop, std_prop)

for (i in c('H3K27ac', 'H3K4me3', 'ATAC', 'DNase', 'H3K27me3', 'H3K4me1', 'H3K9me3', 'H3K36me3')) {
  p <- results %>% filter(assay == i) %>%
    ggplot(aes(x = geneType, y = avg_prop, fill = geneType, col = geneType)) +
    geom_bar(stat = "identity", position = "stack")  +
    labs(x = '', y = 'Proportion of promoters', fill = 'Gene type', color = 'Gene type') +
    theme_minimal() + facet_wrap(~tissue, ncol = 4) +
    geom_errorbar(aes(ymin = avg_prop, ymax = pmin(1,avg_prop + std_prop)), 
                  position = position_dodge(0.6), width = 0.5) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(limits = c('protein_coding', 'lncRNA', 'pseudogene'),
                     labels = c('Protein-coding', 'lncRNA', 'Pseudogene'))+
    scale_color_manual(values = cols_1, labels = c("lncRNA" = "lncRNA", "protein_coding" = "Protein-coding", "pseudogene" = "Pseudogene")) +
    scale_fill_manual(values = cols_1, labels = c("lncRNA" = "lncRNA", "protein_coding" = "Protein-coding", "pseudogene" = "Pseudogene")) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),  
      panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = 'bottom')
  ggsave(paste0('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/tissue-by-tissue/pgene_combined/', i, '.pdf'), p, width = 10, height = 11, units = 'in')
}

## only pseudogenes
# combine three types of pseudogenes, and compare them with protein_coding and lncRNA transcripts
long_data <- df %>% filter(RAMPAGE == "True", str_detect(geneType, 'pseudogene')) %>%
  mutate(tss_length = end - start) %>%
  pivot_longer(cols = c(H3K27ac, H3K27me3, H3K4me1, H3K4me3, H3K9me3, H3K36me3, DNase, ATAC), names_to = 'assay', values_to = 'coverage')

results <- long_data %>% na.omit() %>%
  group_by(geneType, assay, tissue) %>%
  summarize(
    avg_prop = mean(coverage / tss_length),
    std_prop = sd(coverage / tss_length)) %>%
  select(geneType, assay, tissue, avg_prop, std_prop)

for (i in c('H3K27ac', 'H3K4me3', 'ATAC', 'DNase', 'H3K27me3', 'H3K4me1', 'H3K9me3', 'H3K36me3')) {
  p <- results %>% filter(assay == i, geneType != 'unitary_pseudogene') %>%
    ggplot(aes(x = geneType, y = avg_prop, fill = geneType, col = geneType)) +
    geom_bar(stat = "identity", position = "stack")  +
    labs(x = '', y = 'Proportion of promoters', fill = 'Gene type', color = 'Gene type') +
    theme_minimal() + facet_wrap(~tissue, ncol = 4) +
    geom_errorbar(aes(ymin = avg_prop, ymax = pmin(1,avg_prop + std_prop)), 
                  position = position_dodge(0.6), width = 0.5) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(limits = c('unprocessed_pseudogene', 'processed_pseudogene'),
                     labels = c('Unprocessed', 'Processed'))+
    scale_color_manual(values = cols_2, labels = c("unprocessed_pseudogene" = "Unprocessed", "processed_pseudogene" = "Processed")) +
    scale_fill_manual(values = cols_2, labels = c("unprocessed_pseudogene" = "Unprocessed", "processed_pseudogene" = "Processed")) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),  
      panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = 'bottom')
  ggsave(paste0('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/tissue-by-tissue/unprocessed_processed/', i, '.pdf'), p, width = 4, height = 6, units = 'in')
}
