rm(list = ls())
library(dplyr)
library(tidyr)
library(tidyverse)
pgene_mapping <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/parent/tss1kb_comparison/tss1kb_inclusion_identity.txt', sep = '\t', header = TRUE)
epidata <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = TRUE)
pgene_mapping <- pgene_mapping %>% mutate(promoter_inherit = ifelse(Inclusion_Pct > 50 & Identity_Pct > 75, "Inherited", "NotInherited"))
pgene_mapping <- pgene_mapping %>% select(PgeneEnsembl103ID, promoter_inherit) %>% distinct()

pgene_epidata <- epidata %>%
  filter(geneType %in% c("unprocessed_pseudogene", "processed_pseudogene")) %>%
  left_join(pgene_mapping, by = c("transcript_id" = "PgeneEnsembl103ID")) %>% filter(!is.na(promoter_inherit))

long_data <- pgene_epidata %>% filter(RAMPAGE == "True") %>%           
  mutate(tss_length = end - start) %>%
  pivot_longer(cols = c(H3K27ac, H3K27me3, H3K4me1, H3K4me3, H3K9me3, H3K36me3, DNase, ATAC),
    names_to = "assay", values_to = "coverage")

top_tbl <- long_data %>% group_by(gene_id, tissue) %>%
  slice_max(avgTPM, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(gene_id, tissue, transcript_id)

one_iso <- long_data %>%
  semi_join(top_tbl,
            by = c("gene_id", "tissue", "transcript_id"))

df_bar <- one_iso %>%
  filter(
    geneType == 'unprocessed_pseudogene',
    assay %in% c("H3K27ac", "H3K4me3", "DNase", "ATAC")) %>%
  mutate(prop = coverage / tss_length) %>%
  group_by(assay, promoter_inherit) %>%
  summarize(
    mean_prop = mean(prop, na.rm = T),
    sd_prop   = sd(prop, na.rm = T),
    n         = n(),
    .groups = "drop"
  )

write.table(df_bar, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/tables/unprocessed_pgene_inheritance_0.50.txt', quote = F, sep = '\t')

p <- ggplot(df_bar, aes(x = assay, y = mean_prop, fill = promoter_inherit, colour = promoter_inherit)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean_prop, ymax = mean_prop + sd_prop, colour = promoter_inherit), position = position_dodge(width = 0.7), width = 0.2) +
  scale_fill_manual(values = c(Inherited = "#8F2D56", NotInherited = "#218380")) +
  scale_color_manual(values = c(Inherited = "#8F2D56", NotInherited = "#218380")) +
  labs(
    title = 'Inclusion - 0.50',
    x = NULL,
    y = "Proportion of promoters") +
  theme_bw(base_size = 14) + ylim(c(0,1)) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(color="black", size = 12), 
    legend.title = element_blank(),   
    legend.position = c(0.2, 0.9),
    legend.background = element_rect(fill = "white",color = "black", size = 0.3))

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/plots/unprocessed_pgene_prop_histoneMarks_promoter_inherited.pdf', p, width = 6, height = 6)
