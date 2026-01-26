## epigenetic characterization of pseudoegens with and without LINEs/SINEs
rm(list = ls())
library(stringr)
library(dplyr)
library(ggplot2)
library(ComplexUpset)
library(tidyverse)
# overlap with TEs
setwd('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/intersection/repetitive_elements/')
rmsk <- read.table('promoter_regions_rmsk.txt')
rmsk <- rmsk %>% mutate(element_len = if_else(V8 != '.', V10-V9, NA)) %>% 
  dplyr::rename(overlap_len = V23, chr = V1, start = V2, end = V3, gene_id = V4, score = V5, strand = V6, transcript_id = V7)
# promoter table
df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)

plot_histoneMarks <- function(df, name){
  top_tbl <- df %>%
    group_by(gene_id, tissue) %>%                      
    slice_max(avgTPM, n = 1, with_ties = FALSE) %>%    
    ungroup() %>% select(gene_id, tissue, transcript_id)
  one_iso <- df %>% semi_join(top_tbl, by = c("gene_id", "tissue", "transcript_id"))
  results <- one_iso %>% na.omit() %>% filter(assay %in% c('H3K27ac', 'H3K4me3', 'DNase', 'ATAC')) %>% group_by(geneType, assay) %>%
    summarize(
      avg_prop = mean(coverage / tss_length),
      std_prop = sd(coverage / tss_length)) %>% select(geneType, assay, avg_prop, std_prop)
  results$assay <- factor(results$assay, levels = c('H3K27ac', 'H3K4me3', 'ATAC', 'DNase'))
  cols <- c('unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b')
  p <- results %>% filter(geneType != 'unitary_pseudogene') %>%
    ggplot(aes(x=assay, y=avg_prop, fill=geneType)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.8) +
    geom_errorbar(aes(ymin = avg_prop, ymax = pmin(avg_prop + std_prop, 1), color = geneType),
                  position = position_dodge(width = 0.8), width = 0.4, size = 0.5) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    labs(x = "", y = "Proportion of promoters", title = name) +
    theme_bw(base_size = 14) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, size = 1),  
      panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text = element_text(color="black", size = 12),
      legend.position = 'right') 
  
  return(p)
}


pgene_LINE_SINE <- df %>% filter(RAMPAGE== "True") %>%
  left_join(rmsk, by = c('chr', 'start', 'end', 'gene_id', 'score', 'strand','transcript_id')) %>%
  dplyr::rename(repClass = V18) %>%
  filter(str_detect(repClass, 'LINE|SINE'), str_detect(geneType, 'pseudogene')) %>%
  select(chr, start, end, gene_id, transcript_id, geneType, tissue, avgTPM,
         matches("H3K|ATAC|DNase")) %>%
  mutate(tss_length = end - start) %>%
  pivot_longer(cols = c(H3K27ac, H3K27me3, H3K4me1, H3K4me3, H3K9me3, H3K36me3, DNase, ATAC), names_to = 'assay', values_to = 'coverage')

pgene_no_LINE_SINE <- df %>% filter(RAMPAGE== "True") %>%
  left_join(rmsk, by = c('chr', 'start', 'end', 'gene_id', 'score', 'strand','transcript_id')) %>%
  dplyr::rename(repClass = V18) %>%
  filter(!str_detect(repClass, 'LINE|SINE'), str_detect(geneType, 'pseudogene')) %>%
  select(chr, start, end, gene_id, transcript_id, geneType, tissue, avgTPM,
         matches("H3K|ATAC|DNase")) %>%
  mutate(tss_length = end - start) %>%
  pivot_longer(cols = c(H3K27ac, H3K27me3, H3K4me1, H3K4me3, H3K9me3, H3K36me3, DNase, ATAC), names_to = 'assay', values_to = 'coverage')

p1 <- plot_histoneMarks(pgene_LINE_SINE, 'with LINE/SINE')
p2 <- plot_histoneMarks(pgene_no_LINE_SINE, 'without LINE/SINE')

library(patchwork)
p <- (p1 | p2) + plot_layout(guides = 'collect', axes = 'collect')
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/plots/prop_promoters_with_LINEsandSINEs.pdf', p, width = 10, height = 6)
