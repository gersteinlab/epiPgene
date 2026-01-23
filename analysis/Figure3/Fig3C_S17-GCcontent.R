rm(list = ls())
library(stringr)
library(dplyr)
# upstream region
all <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
upstream <- all %>% filter(RAMPAGE ==  "True") %>%
  dplyr::select(chr, start, end, gene_id, score, strand, transcript_id, geneType) %>% distinct() %>%
  mutate(
    new_start = if_else(strand == "-", end - 1000, start),
    new_end = if_else(strand == "-", end, start + 1000),
  ) %>%
  select(-start, -end) %>%
  dplyr::rename(start = new_start, end = new_end) %>%
  dplyr::select(chr, start, end, gene_id, score, strand, transcript_id, geneType)

write.table(upstream, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/GCcontent/GRCh38_tss_upstream.v29.bed', sep = '\t', row.names = F, col.names = F, quote = F)

## downstream
downstream <- all %>% filter(RAMPAGE ==  "True") %>%
  dplyr::select(chr, start, end, gene_id, score, strand, transcript_id, geneType) %>% distinct() %>%
  mutate(
    new_start = if_else(strand == "-", start, start + 1000),
    new_end = if_else(strand == "-", start + 1000, end),
  ) %>%
  select(-start, -end) %>%
  dplyr::rename(start = new_start, end = new_end) %>%
  dplyr::select(chr, start, end, gene_id, score, strand, transcript_id, geneType)

write.table(downstream, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/GCcontent/GRCh38_tss_downstream.v29.bed', sep = '\t', row.names = F, col.names = F, quote = F)

# start from here
library(gghalves)
library(ggpubr)
GC <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/GCcontent/GRCh38_tss_upstream.v29.GC.txt', comment.char = '', header = T)
df <- GC %>% select(X7_usercol, X10_pct_gc) %>% 
  dplyr::rename(transcript_id = X7_usercol, GC_content = X10_pct_gc) %>%
  dplyr::mutate(transcript_id = str_remove(transcript_id, '\\..*'))

cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')
data_to_plot <- all %>% filter(RAMPAGE == "True") %>%
  dplyr::select(chr, start, end, strand, transcript_id, geneType) %>% distinct() %>%
  left_join(df, by = 'transcript_id')
data_to_plot$geneType <- factor(data_to_plot$geneType, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene'))

p <- ggplot(data_to_plot, aes(x = geneType, y = GC_content, fill = geneType))  +
  geom_half_violin(position=position_nudge(x=0.12), side='r', alpha = 0.6, trim = T, aes(color = geneType)) +
  geom_boxplot(width=.15, position=position_nudge(x=0), alpha = 0.6, outliers = F, aes(color = geneType)) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols)+
  stat_compare_means(comparisons = list(c("processed_pseudogene", "protein_coding"), c("processed_pseudogene", "lncRNA"), c("processed_pseudogene", "unprocessed_pseudogene"), c("processed_pseudogene", "unitary_pseudogene")),
                     method = 'wilcox.test', label = 'p.signif', tip.length = 0.01, step.increase = 0.07,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'n.s.'))) +
  labs(x = '', y = 'GC content') +
  scale_x_discrete(labels = c('Protein-coding', 'lncRNA', 'Processed ', 'Unprocessed', 'Unitary'), expand = c(0, 0)) +
  scale_y_continuous(breaks=c(2,3,4,5,6,7,8,9,10)*0.1) +
  theme(axis.text = element_text(colour = 'black'), legend.position = 'None') +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'None',
    axis.text = element_text(color="black", size = 12))
  
tbl <- data_to_plot %>% rstatix::wilcox_test(GC_content~geneType,ref.group = 'processed_pseudogene', p.adjust.method = 'fdr')
write.table(tbl, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/tables/GC_content_upstream_wilcoxon.txt', sep = '\t', quote = F)
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/plots/GC_content_upstream.pdf', p, width = 5, height = 5)

# downstream
GC <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/GCcontent/GRCh38_tss_downstream.v29.GC.txt', comment.char = '', header = T)
df <- GC %>% select(X7_usercol, X10_pct_gc) %>% 
  dplyr::rename(transcript_id = X7_usercol, GC_content = X10_pct_gc) %>%
  dplyr::mutate(transcript_id = str_remove(transcript_id, '\\..*'))

data_to_plot <- all %>% filter(RAMPAGE == "True") %>%
  dplyr::select(chr, start, end, strand, transcript_id, geneType) %>% distinct() %>%
  left_join(df, by = 'transcript_id')
data_to_plot$geneType <- factor(data_to_plot$geneType, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene'))

p <- ggplot(data_to_plot, aes(x = geneType, y = GC_content, fill = geneType))  +
  geom_half_violin(position=position_nudge(x=0.12), side='r', alpha = 0.6, trim = T, aes(color = geneType)) +
  geom_boxplot(width=.15, position=position_nudge(x=0), alpha = 0.6, outliers = F, aes(color = geneType)) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols)+
  stat_compare_means(comparisons = list(c("processed_pseudogene", "protein_coding"), c("processed_pseudogene", "lncRNA"), c("processed_pseudogene", "unprocessed_pseudogene"), c("processed_pseudogene", "unitary_pseudogene")),
                     method = 'wilcox.test', label = 'p.signif', tip.length = 0.01, step.increase = 0.07,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'n.s.'))) +
  labs(x = '', y = 'GC content') +
  scale_x_discrete(labels = c('Protein-coding', 'lncRNA', 'Processed ', 'Unprocessed', 'Unitary'), expand = c(0, 0)) +
  scale_y_continuous(breaks=c(2,3,4,5,6,7,8,9,10)*0.1) +
  theme(axis.text = element_text(colour = 'black'), legend.position = 'None') +
  coord_flip() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'None',
    axis.text = element_text(color="black", size = 12))

tbl <- data_to_plot %>% rstatix::wilcox_test(GC_content~geneType,ref.group = 'processed_pseudogene', p.adjust.method = 'fdr')
write.table(tbl, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/tables/GC_content_downstream_wilcoxon.txt', sep = '\t', quote = F)
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/figS/GC_content_downstream.pdf', p, width = 5, height = 5)