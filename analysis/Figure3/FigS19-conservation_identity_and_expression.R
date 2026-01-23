rm(list = ls())
library(dplyr)
library(rstatix)
library(ggplot2)
cols <- c('Unprocessed' = '#edae49', 'Processed' = '#d1495b')
## phastCons 20-way
phastcons20_conservation <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/phastCons/upstream_phastcons20way.bed')
colnames(phastcons20_conservation) <- c('chr', 'start', 'end', 'id', 'conservation')
mastertable <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/GCcontent/GRCh38_tss_upstream.v29.bed', header = F, sep = '\t')
colnames(mastertable) <- c('chr', 'start', 'end', 'gene_id', 'score', 'strand', 'transcript_id', 'geneType')
phastcons20_conservation <- mastertable %>% left_join(phastcons20_conservation, by = c('chr', 'start', 'end'))

## identity to parent genes
df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/parent/tss1kb_comparison/tss1kb_inclusion_identity.txt', sep = '\t', header = T)

## expression
all <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
expr <- all %>% filter(geneType %in% c('processed_pseudogene', 'unprocessed_pseudogene')) %>% 
  group_by(transcript_id, geneType) %>% summarise(log1p_expr = log1p(median(avgTPM))) %>%
  mutate(geneType = recode(geneType, 'processed_pseudogene' = 'Processed', 'unprocessed_pseudogene' = 'Unprocessed'))


## correlation: phastcons vs. expression
tbl <- phastcons20_conservation %>% merge(expr, by = 'transcript_id') %>%
  group_by(geneType.y) %>%
  summarise(cor = cor.test(conservation, log1p_expr, use = "spearman")$estimate, p_value = cor.test(conservation, log1p_expr, use = "spearman")$p.value)
write.table(tbl, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/tables/promoter_phastcons20_expression_correlation.txt', sep = '\t', quote = F)

p1 <- phastcons20_conservation %>% merge(expr, by = 'transcript_id') %>%
  ggplot(aes(x = log1p_expr, y = conservation)) +
  geom_point(aes(color = geneType.y)) +
  labs(x = 'Transcript expression (logTPM)', y = 'Promoter conservation (phastCons 20-way)') +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(color="black", size = 12), 
    legend.title = element_blank(),   
    legend.position = c(0.85, 0.9),
    legend.background = element_rect(fill = "white",color = "black", size = 0.3))

## parent gene identity
tbl <- df %>% merge(expr, by.x = 'PgeneEnsembl103ID', by.y = 'transcript_id') %>%
  group_by(geneType.y) %>%
  summarise(cor = cor.test(Identity_Pct, log1p_expr, use = "spearman")$estimate, p_value = cor.test(Identity_Pct, log1p_expr, use = "spearman")$p.value)
write.table(tbl, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/tables/promoter_identity_parentGene_expression_correlation.txt', sep = '\t', quote = F)

p2 <- df %>% merge(expr, by.x = 'PgeneEnsembl103ID', by.y = 'transcript_id') %>%
  ggplot(aes(x = log1p_expr, y = Identity_Pct)) +
  geom_point(aes(color = geneType.y)) +
  labs(x = 'Transcript expression (logTPM)', y = 'Promoter sequence identity to parent gene (%)') +
  scale_color_manual(values = cols) +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(color="black", size = 12), 
    legend.title = element_blank(),   
    legend.position = c(0.85, 0.9),
    legend.background = element_rect(fill = "white",color = "black", size = 0.3))

library(ggExtra)
library(patchwork)
library(grid)
library(gridGraphics)

p3 <- ggMarginal(p1, type = 'boxplot', groupColour = T)
p4 <- ggMarginal(p2, type = 'boxplot', groupColour = T)

p5 <- grid.grabExpr(print(p3))
p6 <- grid.grabExpr(print(p4))
library(patchwork)
p <- (wrap_elements(full = p5) + wrap_elements(full = p6)) + plot_layout(axes = 'collect')
p <- (p1 + p2) + plot_layout(axes = 'collect')
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/plots/conservation_identity_and_expression.pdf', p, width = 10, height = 6, units = 'in')
