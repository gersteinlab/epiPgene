rm(list = ls())
library(dplyr)
library(stringr)
library(ggplot2)
library(ggmagnify)
library(tidyverse)
df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
geneType_reordered <- df %>% dplyr::select(tissue, transcript_id, geneType, RAMPAGE) %>% distinct() %>%
  mutate(geneType = if_else(str_detect(geneType, 'pseudogene'), 'pseudogene', geneType)) %>%
  group_by(geneType) %>%
  summarise(RAMPAGE_TRUE = sum(RAMPAGE == "True"), total = n()) %>%
  mutate(proportion = RAMPAGE_TRUE / total) %>%
  arrange(-proportion)

data_to_plot <- df %>% dplyr::select(tissue, transcript_id, geneType, RAMPAGE) %>% distinct() %>% mutate(geneType = if_else(str_detect(geneType, 'pseudogene'), 'pseudogene', geneType))
data_to_plot[['geneType']] <- factor(data_to_plot$geneType, levels = c('protein_coding', 'lncRNA', 'pseudogene'))
data_to_plot[['RAMPAGE']] <- factor(data_to_plot$RAMPAGE, levels = c('False', 'True'))

cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'pseudogene' = '#9392BE')
data_to_plot$FillColor <- as.factor(ifelse(data_to_plot$RAMPAGE == "True", cols[data_to_plot$geneType], "white"))
data_to_plot$FillColor <- factor(data_to_plot$FillColor, levels = c("white", setdiff(levels(data_to_plot$FillColor), "white")))

p1 <- ggplot(data_to_plot, aes(x = geneType, fill =FillColor, col = geneType)) + 
  geom_bar(position = 'fill', width = 0.6) +
  scale_fill_identity() +
  scale_color_manual(values = cols) +
  scale_x_discrete(labels = c("Protein-coding", "lncRNA", "Pseudogene")) +
  labs(x = '', y = 'Proportion overlapping RAMPAGE peaks') + theme_classic() + scale_y_continuous(breaks = seq(0,1,0.2)) +
  geom_text(data = geneType_reordered, aes(x = geneType, y = 1.05, label = total), inherit.aes = F) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'None',
    axis.text = element_text(color="black", size = 12),
    axis.title = element_text(color="black",size = 12))

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig1/fig1C.pdf', p1, width = 6, height = 4)

## stratified by pseudogene biotype
pgeneType_reordered <- df %>% filter(str_detect(geneType, 'pseudogene')) %>%
  dplyr::select(tissue, transcript_id, geneType, RAMPAGE) %>% distinct() %>%
  group_by(geneType) %>%
  summarise(RAMPAGE_TRUE = sum(RAMPAGE == "True"), total = n()) %>%
  mutate(proportion = RAMPAGE_TRUE / total) %>%
  arrange(-proportion)

data_to_plot <- df %>% filter(str_detect(geneType, 'pseudogene')) %>% dplyr::select(tissue, transcript_id, geneType, RAMPAGE) %>% distinct()
data_to_plot[['geneType']] <- factor(data_to_plot$geneType, levels = c('unprocessed_pseudogene', 'processed_pseudogene', 'unitary_pseudogene'))
data_to_plot[['RAMPAGE']] <- factor(data_to_plot$RAMPAGE, levels = c('False', 'True'))
cols <- c('unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')
data_to_plot$FillColor <- as.factor(ifelse(data_to_plot$RAMPAGE == "True", cols[data_to_plot$geneType], "white"))
data_to_plot$FillColor <- factor(data_to_plot$FillColor, levels = c("white", setdiff(levels(data_to_plot$FillColor), "white")))

p2 <- ggplot(data_to_plot, aes(x = geneType, fill =FillColor, col = geneType)) + 
  geom_bar(position = 'fill', width = 0.6) +
  scale_fill_identity() +
  scale_color_manual(values = cols) +
  scale_x_discrete(labels = c("Unprocessed", "Processed", "Unitary")) +
  labs(x = '', y = 'Proportion overlapping RAMPAGE peaks') + theme_classic() + scale_y_continuous(breaks = seq(0,1,0.2)) +
  geom_text(data = pgeneType_reordered, aes(x = geneType, y = 1.05, label = total), inherit.aes = F) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'None',
    axis.text = element_text(color="black", size = 12),
    axis.title = element_text(color="black",size = 12)) + coord_flip()

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig1/fig1D.pdf', p2, width = 6, height = 4)
#
write.table(geneType_reordered, '/gpfs/gibbs/pi/gerstein/yj329/epi/analysis/Fig1/rPeaks_all.txt', quote = F, row.names = F)
write.table(pgeneType_reordered, '/gpfs/gibbs/pi/gerstein/yj329/epi/analysis/Fig1/rPeaks_pgene.txt', quote = F, row.names = F)

# stats
matrix(c(100667, 134022-100667, 1022, 2879-1022), byrow = T, nrow = 2) %>% as.table() %>% fisher.test() # pcg and lncRNA
matrix(c(100667, 134022-100667, 1664, 4457-1664), byrow = T, nrow = 2) %>% as.table() %>% fisher.test() # pcg and pseudogene
matrix(c(1028, 2680-1028, 594, 1690-594), byrow = T, nrow = 2) %>% as.table() %>% fisher.test() # pcg and pseudogene