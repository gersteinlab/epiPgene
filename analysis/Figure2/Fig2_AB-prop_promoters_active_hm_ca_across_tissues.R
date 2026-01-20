rm(list = ls())
library(dplyr)
library(stringr)
library(reshape2)
library(tidyverse)
df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
cols_1 <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'pseudogene' = '#BEB6DD')
cols_2 <- c('protein_coding' = 'black', 'lncRNA' = 'gray','unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')

# combine three types of pseudogenes, and compare them with protein_coding and lncRNA transcripts
long_data <- df %>% filter(RAMPAGE == "True") %>%
  mutate(geneType = if_else(str_detect(geneType, 'pseudogene'), 'pseudogene', geneType)) %>%
  mutate(tss_length = end - start) %>%
  pivot_longer(cols = c(H3K27ac, H3K27me3, H3K4me1, H3K4me3, H3K9me3, H3K36me3, DNase, ATAC), names_to = 'assay', values_to = 'coverage')

top_tbl <- long_data %>%                              
  group_by(gene_id, tissue) %>%                       
  slice_max(avgTPM, n = 1, with_ties = FALSE) %>%     
  ungroup() %>%                                       
  select(gene_id, tissue, transcript_id)

one_iso <- long_data %>% 
  semi_join(top_tbl, by = c("gene_id", "tissue", "transcript_id"))

# proportion across tissues
results <- one_iso %>% na.omit() %>% filter(assay %in% c('H3K27ac', 'H3K4me3', 'DNase', 'ATAC')) %>%
  group_by(geneType, assay) %>%
  summarize(
   avg_prop = mean(coverage / tss_length),
   std_prop = sd(coverage / tss_length)) %>%
  select(geneType, assay, avg_prop, std_prop)

results$assay <- factor(results$assay, levels = c('H3K27ac', 'H3K4me3', 'ATAC', 'DNase'))
p <- ggplot(results, aes(x = geneType, y = avg_prop, fill = geneType, col = geneType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = '', y = 'Proportion') +
  theme_minimal() + facet_wrap(~ assay, nrow = 1) +
  geom_errorbar(aes(ymin = avg_prop, ymax = pmin(1,avg_prop + std_prop)), 
                position = position_dodge(0.6), width = 0.5) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(limits = c('protein_coding', 'lncRNA', 'pseudogene'),
    labels = c('Protein-coding', 'lncRNA', 'Pseudogene'))+
  scale_color_manual(values = cols_1) +
  scale_fill_manual(values = cols_1) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text = element_text(color="black", size = 12),
    legend.position = 'None')
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig2/proportion_across_tissue.pdf', p, width = 8, height = 4)

# only pseudogenes
# combine three types of pseudogenes, and compare them with protein_coding and lncRNA transcripts
long_data <- df %>% filter(RAMPAGE == 'True', str_detect(geneType, 'pseudogene')) %>%
  mutate(tss_length = end - start) %>%
  pivot_longer(cols = c(H3K27ac, H3K27me3, H3K4me1, H3K4me3, H3K9me3, H3K36me3, DNase, ATAC), names_to = 'assay', values_to = 'coverage')

top_tbl <- long_data %>%                              
  group_by(gene_id, tissue) %>%                       
  slice_max(avgTPM, n = 1, with_ties = FALSE) %>%     
  ungroup() %>%                                       
  dplyr::select(gene_id, tissue, transcript_id)

one_iso <- long_data %>% 
  semi_join(top_tbl, by = c("gene_id", "tissue", "transcript_id"))

results <- one_iso %>% na.omit() %>% filter(assay %in% c('H3K27ac', 'H3K4me3', 'DNase', 'ATAC')) %>%
  group_by(geneType, assay) %>%
  summarize(
    avg_prop = mean(coverage / tss_length),
    std_prop = sd(coverage / tss_length)) %>%
  dplyr::select(geneType, assay, avg_prop, std_prop)

results$assay <- factor(results$assay, levels = c('H3K27ac', 'H3K4me3', 'ATAC', 'DNase'))
p <- ggplot(results, aes(x = geneType, y = avg_prop, fill = geneType, col = geneType)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_errorbar(aes(ymin = avg_prop, ymax = pmin(1,avg_prop + std_prop)), 
                position = position_dodge(0.6), width = 0.5) +
  labs(x = '', y = 'Proportion') +
  theme_minimal() + facet_wrap(~ assay, nrow = 1) + 
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(limits = c('unitary_pseudogene', 'unprocessed_pseudogene', 'processed_pseudogene'),
                   labels = c('Unitary', 'Unprocessed', 'Processed'))+
  scale_color_manual(values = cols_2) +
  scale_fill_manual(values = cols_2) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1), 
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
    axis.text = element_text(color="black", size = 12),
    legend.position = 'None')
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig2/proportion_across_tissue_pgene.pdf', p, width = 8, height = 4)