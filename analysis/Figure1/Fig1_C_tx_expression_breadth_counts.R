rm(list = ls())
setwd('/gpfs/gibbs/pi/gerstein/yj329/epiPgene')

df <- read.table('freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', header = T, sep = '\t')
cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'pseudogene' = '#9392BE')

tbl <- df %>%
  mutate(geneType = if_else(str_detect(geneType, 'pseudogene'), 'pseudogene', geneType)) %>%
  group_by(geneType, transcript_id) %>%
  summarise(num_tissues = n_distinct(tissue)) %>%
  mutate(
  num_tissues_group = cut(
    num_tissues,
    breaks = c(-Inf, 1, 2, 3, 4, 5, Inf),
    labels = c("1", "2", "3", "4", "5", "6+")
  )) %>%
  select(-num_tissues) %>%
  group_by(geneType, num_tissues_group) %>% 
  count(num_tissues_group, name = "num_transcripts")

tbl$geneType <- factor(tbl$geneType, levels = c('pseudogene', 'protein_coding', 'lncRNA'))

# Plot
ggplot(tbl, aes(x = num_tissues_group, y = num_transcripts, fill = geneType)) +
  geom_bar(stat = "identity") +
  facet_wrap(~geneType, nrow = 3, scales = 'free_y', labeller = as_labeller(function(x) {
    x <- gsub("pseudogene", "Pseudogene", x)
    x <- gsub("protein_coding", "Protein-coding", x)
    return(x)
  })) + 
  scale_fill_manual(values = cols) + 
  labs(x = "Number of tissues", y = "Number of expressed transcripts") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'None',
    axis.text = element_text(color="black", size = 12),
    axis.title = element_text(color="black",size = 12))

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig1/fig1B.pdf', width = 4, height = 6)
