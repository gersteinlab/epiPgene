setwd('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/HiC_interval/')
rm(list = ls())
library(dplyr)
library(tidyr)
library(dplyr)
ccre <- read.table('hic_interveals_w_cCREs.txt')
ccre_processed <- ccre %>% group_by(V1, V2, V3, V9) %>%
  summarise(ccre_len = sum(V10)) %>% ungroup() %>%
  complete(nesting(V1, V2, V3), V9, fill = list(ccre_len = 0)) %>% filter(V9 != '.') %>% mutate(tot_len = V3 - V2) %>%
  ungroup() %>%
  pivot_wider(names_from = V9, values_from = ccre_len) %>% 
  mutate(across(c(PLS, TF, dELS, pELS, CA, `CA-CTCF`, `CA-H3K4me3`, `CA-TF`), ~./tot_len, .names = "{col}"))

df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig4/tables/hic_coords_transcripts.txt')

data_to_plot <- df %>% dplyr::select(V1, V2, V3, V4, V5, V6, V7, V8, V13) %>% 
  merge(.,y = ccre_processed, by.x = c('V6', 'V7', 'V8'), by.y = c('V1', 'V2', 'V3')) %>%
  dplyr::rename(hic_chr = V6, hic_start = V7, hic_end = V8) %>% 
  group_by(hic_chr, hic_start, hic_end) %>%
  filter(!duplicated(V13)) %>% ## remove if a given HiC interval always intersect the same gene type
  ungroup()

# chromatin accessibility
library(stringr)
tmp <- data_to_plot %>% pivot_longer(cols = c(PLS, TF, dELS, pELS, CA, `CA-CTCF`, `CA-H3K4me3`, `CA-TF`), names_to = "attr", values_to = "Value") %>%
  filter(str_detect(attr, 'CA')) %>%
  group_by(V13, attr) %>%
  summarise(Mean = mean(Value, na.rm = TRUE),
            SE = sd(Value, na.rm = TRUE)/sqrt(n()))

cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b')
library(ggplot2)
tmp$V13 <- factor(tmp$V13, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene'))
p <- ggplot(tmp, aes(x=attr, y=Mean, fill=V13)) + 
  geom_bar(stat="identity", position=position_dodge(), alpha = 1) +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE, color = V13), width = .5, position=position_dodge(.9), size = 0.5) +
  scale_fill_manual(values = cols, labels = c("protein_coding" = "Protein-coding", "lncRNA" = "lncRNA", "processed_pseudogene" = "Processed", "unprocessed_pseudogene" = "Unprocessed")) +
  scale_color_manual(values = cols, labels = c("protein_coding" = "Protein-coding", "lncRNA" = "lncRNA", "processed_pseudogene" = "Processed", "unprocessed_pseudogene" = "Unprocessed")) +
  labs(x="", y="Proportion", fill = 'Biotype', color = 'Biotype') +
  theme(
    panel.background = element_rect(fill = "white", color = NA),  
    plot.background = element_rect(fill = "white", color = NA),   
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(color="black", size = 12),
    legend.position = 'right')

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/figS/hic_interval_cCRE_CA.pdf', p , width = 8, height = 6)

# ELS, etc
tmp <- data_to_plot %>% pivot_longer(cols = c(PLS, TF, dELS, pELS, CA, `CA-CTCF`, `CA-H3K4me3`, `CA-TF`), names_to = "attr", values_to = "Value") %>%
  filter(!str_detect(attr, 'CA')) %>%
  group_by(V13, attr) %>%
  summarise(Mean = mean(Value, na.rm = TRUE),
            SE = sd(Value, na.rm = TRUE)/sqrt(n()))
cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b')
tmp$V13 <- factor(tmp$V13, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene'))
p <- ggplot(tmp, aes(x=attr, y=Mean, fill=V13)) + 
  geom_bar(stat="identity", position=position_dodge(), alpha = 1) +
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE, color = V13), width = .5, position=position_dodge(.9), size = 0.5) +
  scale_fill_manual(values = cols, labels = c("protein_coding" = "Protein-coding", "lncRNA" = "lncRNA", "processed_pseudogene" = "Processed", "unprocessed_pseudogene" = "Unprocessed")) +
  scale_color_manual(values = cols, labels = c("protein_coding" = "Protein-coding", "lncRNA" = "lncRNA", "processed_pseudogene" = "Processed", "unprocessed_pseudogene" = "Unprocessed")) +
  labs(x="", y="Proportion",fill = 'Biotype', color = 'Biotype') +
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA),  
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(color="black", size = 12),
    legend.position = 'right')

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig4/plots/hic_interval_cCRE_others.pdf', p, width = 8, height = 6)

## TSS
tss <- read.table('hic_interveals_w_TSSs.txt')
tss_procesed <- tss %>% group_by(V1, V2, V3) %>%
  mutate(is_overlapping_any_tss = if_else(V4 != '.', 1, 0)) %>% ungroup() %>%
  dplyr::select(V1, V2, V3, is_overlapping_any_tss) %>% distinct()

data_to_plot <- df %>% dplyr::select(V1, V2, V3, V4, V5, V6, V7, V8, V13) %>% 
  merge(.,y = tss_procesed, by.x = c('V6', 'V7', 'V8'), by.y = c('V1', 'V2', 'V3')) %>%
  dplyr::rename(hic_chr = V6, hic_start = V7, hic_end = V8) %>% 
  group_by(hic_chr, hic_start, hic_end) %>%
  filter(!duplicated(V13)) %>% ## remove if a given HiC interval always intersect the same gene type
  ungroup()

tmp <- data_to_plot %>% pivot_longer(cols = is_overlapping_any_tss, names_to = "attr", values_to = "Value")
tmp$V13 <- factor(tmp$V13, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene'))
library(ggmosaic)
library(ggsci)
p <- tmp %>% mutate(V13 = recode(V13, "protein_coding" = "Protein-coding", "lncRNA" = "lncRNA", "processed_pseudogene" = "Processed", "unprocessed_pseudogene" = "Unprocessed")) %>%
  ggplot() +
  geom_mosaic(aes(weight = 1, x = product(V13), fill = Value), show.legend = FALSE) +
  scale_fill_tron() +
  theme_mosaic() +
  facet_grid(~V13) +
  theme(aspect.ratio = 3,
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(x = '', y = 'Overlap any TSS')

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/figS/hic_interval_overlapping_TSS.pdf', p , width = 10, height = 5)
count <- table(tmp$V13, tmp$Value)
write.table(count, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig4/tables/hic_interval_overlapping_TSSs.txt', sep = '\t', quote = F)