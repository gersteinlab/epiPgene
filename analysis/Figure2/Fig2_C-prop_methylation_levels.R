rm(list = ls())
library(stringr)
library(dplyr)
library(ggplot2)
setwd('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/intersection/DNAme/')
files <- list.files(path = '.', pattern = "\\.txt$", full.names = TRUE)

DNAme <- data.frame()
for(file in files) {
  tissueName <- str_split(basename(file), '\\.')[[1]][1]
  temp <- read.table(file, header = TRUE, sep = "\t") %>% mutate(tissue = tissueName)
  DNAme <- rbind(DNAme, temp)
}
DNAme$transcript_id <- gsub("\\..*", "", DNAme$transcript_id)

df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
df2 <- df %>% left_join(DNAme, by = c('transcript_id', 'tissue'))
cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')

## number of CpG sites
median_num_CpG <- df2 %>% dplyr::select(transcript_id, geneType, num_of_CpGs) %>% distinct() %>% pull(num_of_CpGs) %>% median()
p <- df2 %>% dplyr::select(transcript_id, geneType, num_of_CpGs) %>% distinct() %>%
  ggplot(aes(x = geneType, y = as.numeric(num_of_CpGs))) + geom_boxplot(aes(fill = geneType, color = geneType), alpha = 0.6) +
  scale_x_discrete(limits = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene'), 
                   labels = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed', 'Unitary')) +
  scale_fill_manual(values = cols, name = '') +
  scale_color_manual(values = cols, name = '') +
  geom_hline(yintercept = median_num_CpG, linetype = "dashed", color = "#023047") +
  labs(x = '', y = 'Number of CpG sites') + 
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'none',
    axis.text = element_text(color="black", size = 12))
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/intersection_up1kb_only/num_CpG_islands.pdf', width = 6, height = 5)

## methylation level
p <- df2 %>% filter(RAMPAGE == "True") %>%
  dplyr::select(transcript_id, geneType, tissue, num_of_CpGs, Heterogeneously.methylated, Hypermethylated, Hypomethylated) %>% distinct() %>%
  filter(num_of_CpGs != 0) %>%
  group_by(geneType) %>%
  summarise(
    tot_heterogeneous = sum(Heterogeneously.methylated),
    tot_hyper = sum(Hypermethylated),
    tot_hypo = sum(Hypomethylated),
    N = sum(num_of_CpGs)
  ) %>%
  mutate(
    prop_heterogeneous = tot_heterogeneous / N,
    prop_hyper = tot_hyper / N,
    prop_hypo = tot_hypo / N
  ) %>%
  reshape2::melt(id.vars = "geneType", 
                  measure.vars = c("prop_heterogeneous", "prop_hyper", "prop_hypo")) %>%
  ggplot(aes(x=geneType, y = value, 
             fill = factor(variable, levels = c('prop_hypo', 'prop_heterogeneous', 'prop_hyper')))) +
  geom_bar(stat = "identity", color = 'black', size = 0.3) + 
  scale_fill_manual(
    name = 'CpG methylation',
    values = c("prop_hypo" = "#a3b18a", "prop_heterogeneous" = "#588157","prop_hyper" = "#3a5a40"),
    labels = c("Hypomethylated", "Heterogeneously methylated", "Hypermethylated")) +
  labs(y = 'Proportion', x = '') +
  scale_x_discrete(limits = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene'), 
                   labels = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed', 'Unitary')) +
  theme_classic()  +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'bottom',
    axis.text = element_text(color="black", size = 12))
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig2/prop_CpG_methylation.pdf', p, width = 7, height = 4)