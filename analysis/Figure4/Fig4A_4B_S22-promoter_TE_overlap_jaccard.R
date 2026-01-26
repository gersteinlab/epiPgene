rm(list = ls())
library(stringr)
library(dplyr)
library(ggplot2)
library(ComplexUpset)
library(tidyverse)
setwd('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/intersection/repetitive_elements/')
rmsk <- read.table('promoter_regions_rmsk.txt')
rmsk <- rmsk %>% mutate(element_len = if_else(V8 != '.', V10-V9, NA)) %>% 
  dplyr::rename(overlap_len = V23, chr = V1, start = V2, end = V3, gene_id = V4, score = V5, strand = V6, transcript_id = V7)

mapping <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/fig3A_tss1kb_comparison/tss1kb_inclusion_identity.txt', sep = '\t', header = T)
parentGene_unprocessed <- mapping %>% filter(geneType == 'Unprocessed') %>% pull(ParentGene) %>% unique()
parentGene_processed <- mapping %>% filter(geneType == 'Processed') %>% pull(ParentGene) %>% unique()

df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
res <- df %>% filter(RAMPAGE== "True") %>%
  select(chr, start, end, gene_id, score, strand, transcript_id, geneType) %>%
  distinct() %>% left_join(rmsk, by = c('chr', 'start', 'end', 'gene_id', 'score', 'strand','transcript_id')) %>%
  dplyr::rename(repClass = V18)

data_to_plot <- res %>% 
  filter(!str_detect(repClass, "\\?")) %>%
  filter(repClass != 'Unknown') %>%
  mutate(repClass_simplified = str_replace(repClass, "/.*", ""), repClass_simplified = if_else(repClass_simplified == ".", "NA", repClass_simplified)) %>%
  select(transcript_id, geneType, repClass_simplified) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = repClass_simplified, values_from = value, values_fill = list(value = 0)) %>%
  dplyr::rename('Low complexity' = Low_complexity, 'Simple repeat' = Simple_repeat) %>% as.data.frame() %>% select(-'NA')


d1 <- data_to_plot %>% filter(geneType != 'protein_coding')  
d2 <- data_to_plot %>% filter(geneType == 'protein_coding')  

d1 <- d1 %>% mutate(geneType = case_when(
  geneType == "processed_pseudogene" ~ "Processed",
  geneType == "unprocessed_pseudogene" ~"Unprocessed",
  geneType == "unitary_pseudogene" ~ "Unitary",
  TRUE ~ geneType
))
  
d1 %>%
  group_by(geneType) %>%
  summarise(
    LINE_only = sum(LINE == 1 & SINE == 0),
    SINE_only = sum(SINE == 1 & LINE == 0),
    LINE_SINE_both = sum(LINE == 1 & SINE == 1)
  )

res %>% select(transcript_id, geneType) %>% distinct() %>%  group_by(geneType) %>% summarise(n=n())

p1 <- ComplexUpset::upset(d1, intersect = colnames(data_to_plot)[-c(1,2)], 
                    base_annotations=list(
                      'Intersection size' = intersection_size(
                        counts = FALSE, mapping = aes(fill = factor(geneType, levels = c('lncRNA', 'Processed', 'Unprocessed', 'Unitary')))) + 
                        scale_fill_manual(values = c('lncRNA' = 'gray', 'Unprocessed' = '#edae49', 'Processed' = '#d1495b', 'Unitary' = '#00798c')) +
                        labs(fill = 'Biotype')), 
                    min_size = 10, width_ratio = 0.1,
                    set_sizes = (upset_set_size(
                    geom=geom_bar(fill = 'black', width=0.8), position='right')), guides = 'over')
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig4/plots/upset_nonPCG_repetitive_elements.pdf',p1, width = 12, height = 6)

p2 <- ComplexUpset::upset(d2, intersect = colnames(data_to_plot)[-c(1,2)], 
                          base_annotations=list(
                            'Intersection size' = intersection_size(
                              counts = FALSE, mapping = aes(fill = factor(geneType))) + 
                              scale_fill_manual(values = c('protein_coding' = 'black')) +
                              labs(fill = 'Biotype')), 
                          min_size = 50, width_ratio = 0.1,
                          set_sizes = (upset_set_size(
                          geom=geom_bar(fill = 'black', width=0.8), position='right')), guides = 'over')
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig4/upset_PCG_repetitive_elements.pdf',p2, width = 12, height = 6)

# jaccard
data_to_plot <- res %>% 
  filter(V8!= '.', !str_detect(repClass, "\\?"),repClass != 'Unknown') %>%
  mutate(repClass_simplified = str_replace(repClass, "/.*", "")) %>%
  select(transcript_id, geneType, overlap_len, element_len, repClass_simplified) %>%
  group_by(transcript_id, geneType, repClass_simplified) %>%
  summarise(overlap = sum(overlap_len), element = sum(element_len)) %>%
  mutate(js = overlap /(2000 + element - overlap))

data_to_plot$geneType <- factor(data_to_plot$geneType)
## jaccard score
cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b')
tmp <- data_to_plot %>% na.omit() %>% filter(geneType != 'unitary_pseudogene', repClass_simplified == 'SINE' | repClass_simplified == 'LINE')

ggplot(tmp, aes(x = geneType, y = js)) +
  geom_violin(aes(fill = geneType), size = 0.3, alpha = 0.6, trim = F) +
  geom_boxplot(width = 0.1, outlier.shape = NA, size = 0.3) +
  scale_fill_manual(values = cols) +
  facet_wrap(~ repClass_simplified, ncol = 1) +
  theme_bw() + labs(x = '', y = 'Jaccard index') +
  scale_x_discrete(limits = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene'),
                   labels = c('Protein-coding', 'lncRNA', 'Processed','Unprocessed')) +
  stat_compare_means(comparisons = list(c("processed_pseudogene", "protein_coding"), c("processed_pseudogene", "lncRNA"), c("processed_pseudogene", "unprocessed_pseudogene")),
                     method = 'wilcox.test', label = 'p.signif', tip.length = 0.01, step.increase = 0.1,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'n.s.'))) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'None',
    axis.text = element_text(color="black", size = 12))
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig4/plots/jaccard.pdf', width = 6, height = 6)
  
library(rstatix)
tbl1 <- data_to_plot %>%
  na.omit() %>%
  filter(geneType != 'unitary_pseudogene', repClass_simplified %in% c('LINE', 'SINE', 'Simple_repeat', 'Low_complexity', 'LTR', 'DNA')) %>% 
  group_by(repClass_simplified) %>% 
  kruskal_test(js ~ geneType) %>%
  adjust_pvalue(method = 'fdr')

tbl2 <- data_to_plot %>%
  na.omit() %>%
  filter(geneType %in% c('processed_pseudogene', 'unprocessed_pseudogene', 'protein_coding', 'lncRNA'), repClass_simplified %in% c('LINE', 'SINE')) %>%
  group_by(repClass_simplified) %>%
  wilcox_test(js~geneType, ref.group = 'processed_pseudogene', p.adjust.method = 'fdr') %>% ungroup()

write.table(tbl1, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig4/tables/jaccard_kw.txt', sep = '\t', quote = F, row.names = F)
write.table(tbl2, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig4/tables/jaccard_post_hoc.txt', sep = '\t', quote = F, row.names = F)