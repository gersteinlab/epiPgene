library(ggplot2)
library(dplyr)
library(tidyverse)
df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
mapping <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/parent/tss1kb_comparison/tss1kb_inclusion_identity.txt', sep = '\t', header = T)
parentGene_unprocessed <- mapping %>% filter(geneType == 'Unprocessed') %>% pull(ParentGene) %>% unique()
parentGene_processed <- mapping %>% filter(geneType == 'Processed') %>% pull(ParentGene) %>% unique()
protein_coding_tx <- df %>% filter(RAMPAGE== "True") %>% filter(geneType == 'protein_coding') %>% pull(transcript_id) %>% unique()
processed_tx <- df %>% filter(RAMPAGE== "True") %>% filter(geneType == 'processed_pseudogene') %>% pull(transcript_id) %>% unique()
unprocessed_tx <- df %>% filter(RAMPAGE== "True") %>% filter(geneType == 'unprocessed_pseudogene') %>% pull(transcript_id) %>% unique()
# parent
parentGene_processed <- mapping %>% filter(PgeneEnsembl103ID %in% processed_tx) %>% pull(ParentGene) %>% unique()
parent_processed_tx <-  df %>% filter(gene_id %in% parentGene_processed) %>% pull(transcript_id) %>% unique()
parentGene_unprocessed <- mapping %>% filter(PgeneEnsembl103ID %in% unprocessed_tx) %>% pull(ParentGene) %>% unique()
parent_unprocessed_tx <-  df %>% filter(gene_id %in% parentGene_unprocessed) %>% pull(transcript_id) %>% unique()
parent_all_tx <- union(parent_processed_tx, parent_unprocessed_tx) %>% unique()
# overlap
setwd('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/intersection/repetitive_elements/')
rmsk <- read.table('promoter_regions_rmsk.txt')
rmsk <- rmsk %>% mutate(element_len = if_else(V8 != '.', V10-V9, NA)) %>% 
  dplyr::rename(overlap_len = V23, chr = V1, start = V2, end = V3, gene_id = V4, score = V5, strand = V6, transcript_id = V7)

res <- df %>%
  select(chr, start, end, gene_id, score, strand, transcript_id, geneType) %>%
  distinct() %>% left_join(rmsk, by = c('chr', 'start', 'end', 'gene_id', 'score', 'strand','transcript_id')) %>%
  dplyr::rename(repClass = V18)

data_to_plot <- res %>% 
  filter(!str_detect(repClass, "\\?")) %>%
  filter(repClass != "Unknown") %>%
  mutate(
    repClass_simplified = str_replace(repClass, "/.*", ""),
    repClass_simplified = if_else(repClass_simplified == ".", "NA", repClass_simplified)
  ) %>%
  select(transcript_id, repClass_simplified) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = repClass_simplified, values_from = value, values_fill = list(value = 0)) %>%
  rename(`Low complexity` = Low_complexity, `Simple repeat` = Simple_repeat) %>%
  as.data.frame() %>%
  select(-`NA`)


pc_minus_parents_tx <- setdiff(protein_coding_tx, parent_all_tx)

group_map <- bind_rows(
  tibble(transcript_id = protein_coding_tx, group = "All_protein_coding"),
  tibble(transcript_id = pc_minus_parents_tx, group = "Protein_coding_minus_parents"),
  tibble(transcript_id = parent_all_tx, group = "Parent_protein_coding_genes"),
  tibble(transcript_id = parent_processed_tx, group = "Parent_for_processed"),
  tibble(transcript_id = parent_unprocessed_tx, group = "Parent_for_unprocessed"),
  tibble(transcript_id = processed_tx, group = "Processed"),
  tibble(transcript_id = unprocessed_tx, group = "Unprocessed")
) %>% distinct()


dat_onehot <- group_map %>%
  left_join(data_to_plot, by = "transcript_id")

plot_df <- dat_onehot %>%
  pivot_longer(
    cols = -c(transcript_id, group),
    names_to = "TE_class",
    values_to = "present"
  ) %>% filter(TE_class %in% c("SINE","LINE","DNA","Simple repeat","LTR"), group != 'Protein_coding_minus_parents') %>%
  group_by(group, TE_class) %>%
  summarise(prop = mean(present, na.rm = TRUE), n = n(), .groups = "drop") %>%
  mutate(
    comparison = dplyr::case_when(
      group %in% c("Processed","Parent_for_processed") ~ "Processed vs parent",
      group %in% c("Unprocessed","Parent_for_unprocessed") ~ "Unprocessed vs parent",
      group %in% c("All_protein_coding","Parent_protein_coding_genes") ~ "All protein-coding vs parents"
    ), group = factor(group, levels = c(
      "Processed","Parent_for_processed",
      "Unprocessed","Parent_for_unprocessed",
      "All_protein_coding","Parent_protein_coding_genes"
    )))

plot_df$TE_class <- factor(plot_df$TE_class, levels = c('LINE', 'SINE', 'LTR', 'DNA', 'Simple repeat'))

library(ggsci)
p <- ggplot(plot_df, aes(x = group, y = prop, fill = TE_class)) +
  geom_col(position = "dodge") +
  coord_flip() +
  facet_wrap(~comparison, ncol = 3, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "% promoters overlapping TE", fill = "TE class") +
  theme_bw(base_size = 11) + scale_fill_manual(values = c('DNA' = '#846267', 'LTR' = '#aeb4a9', 'LINE' = '#00a8e8', 'SINE' = '#90a955', 'Simple repeat' = '#e0c1b3'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color="black", size = 12))
 ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig4/plots/TEoverlap_pgene_parent.pdf', p, width = 14, height = 3)
 
 ## prop test
 ## Processed - LINE
 prop.test(x = c(0.3933333*151, 0.125*72), n = c(151, 72), correct = FALSE) # p = 4.8e-5
 ## Processed - SINE
 prop.test(x = c(0.4*151, 0.472222222*72), n = c(151, 72), correct = FALSE) # p = 0.30
 
 ## Unprocessed - LINE
 prop.test(x = c(0.209580838*167, 0.261904762*42), n = c(167, 42), correct = FALSE) # p = 0.46
 ## Unprocessed - SINE
 prop.test(x = c(0.502994012*167, 0.380952381*42), n = c(167, 42), correct = FALSE) # p = 0.16