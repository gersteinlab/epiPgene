rm(list = ls())
library(patchwork)
setwd('/gpfs/gibbs/pi/gerstein/yj329/epiPgene')
df <- read.table('freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
cols <- c('protein_coding' = 'black', 'lncRNA'='gray' ,'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')

p1 <- df %>%
  filter(geneType == 'protein_coding') %>%
  group_by(geneType, tissue) %>%
  summarise(transcript_count = n()) %>%
  ggplot(aes(x = tissue, y = transcript_count, fill = geneType)) +
  geom_bar(stat = "identity") +
  labs(x = '', y = '# of expressed transcripts', fill = 'geneType') +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1))

p2 <- df %>%
  filter(geneType == 'lncRNA') %>%
  group_by(geneType, tissue) %>%
  summarise(transcript_count = n()) %>%
  # arrange(desc(transcript_count)) %>%
  ggplot(aes(x = tissue, y = transcript_count, fill = geneType)) +
  geom_bar(stat = "identity") +
  labs(x = '', y = '# of expressed transcripts', fill = 'geneType') +
  scale_fill_manual(values = cols) +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 12),
        axis.text.x = element_text(angle = 90, hjust = 1))

data_to_plot <- df %>%
  filter(str_detect(geneType, 'pseudogene')) %>%
  group_by(tissue, geneType) %>%
  summarise(transcript_count = n())

#tmp <- data_to_plot %>%
#  group_by(tissue) %>%
#  summarise(total_count = sum(transcript_count)) %>%
#  arrange(desc(total_count))

p3 <- data_to_plot %>%
  ggplot(aes(x = tissue, y = transcript_count, fill = factor(geneType, levels = c('processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene')))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols) +
  labs(x = '', y = '# of expressed transcripts', fill = 'geneType') +
  theme_bw() +
  theme(axis.text = element_text(colour = 'black', size = 12), axis.text.x = element_text(angle = 90, hjust = 1))

p <- (p1 / p2  / p3) + plot_layout(axes  = "collect", guides = 'collect') & theme(legend.position = 'bottom')
ggsave('num_tx_across_tissue.pdf', p, width = 10, height = 8)

