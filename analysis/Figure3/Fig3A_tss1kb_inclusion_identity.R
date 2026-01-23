library(ggplot2)
library(ggExtra)
library(rstatix)
library(scales)
df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/parent/tss1kb_comparison/tss1kb_inclusion_identity.txt', sep = '\t', header = T)
all <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
df2 <- df %>% mutate(is_subset = PgeneEnsembl103ID %in% all$transcript_id)
cols <- c('Unprocessed' = '#edae49', 'Processed' = '#d1495b')

# expressed
p <- df2 %>% filter(is_subset == TRUE) %>% ggplot(aes(x = Inclusion_Pct, y = Identity_Pct)) +
  geom_point(alpha = 0.8, size = 2, aes(color = geneType)) + 
  scale_x_continuous(trans = pseudo_log_trans(sigma = 20), limits = c(0, 100)) +
  scale_y_continuous(limits = c(60, 100)) +
  labs(x = 'Inclusion (%)', y = 'Identity (%)', main = 'Expressed') +
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
    legend.position = c(0.8, 0.1),
    legend.background = element_rect(
      fill = "white",color = "black", size = 0.3))

p_marg <- ggMarginal(p, type = "boxplot", groupColour = T)
pdf("/gpfs/gibbs/pi/gerstein/yj329/epiPgene/parent/tss1kb_comparison/inclusion_identity_expressed.pdf", width = 6, height = 6)
print(p_marg)
dev.off()

df2 %>% filter(is_subset == TRUE) %>% wilcox_test(Inclusion_Pct~geneType)
df2 %>% filter(is_subset == TRUE) %>% wilcox_test(Identity_Pct~geneType)

# not expressed
p <- df2 %>% filter(is_subset == FALSE) %>% ggplot(aes(x = Inclusion_Pct, y = Identity_Pct)) +
  geom_point(alpha = 0.8, size = 2, aes(color = geneType)) + 
  scale_x_continuous(trans = pseudo_log_trans(sigma = 20), limits = c(0, 100)) +
  scale_y_continuous(limits = c(60, 100)) +
  labs(x = 'Inclusion (%)', y = 'Identity (%)', main = 'Expressed') +
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
    legend.position = c(0.8, 0.1),
    legend.background = element_rect(
      fill = "white",color = "black", size = 0.3))

p_marg <- ggMarginal(p, type = "boxplot", groupColour = T)

pdf("/gpfs/gibbs/pi/gerstein/yj329/epiPgene/parent/tss1kb_comparison/inclusion_identity_notExpressed.pdf", width = 6, height = 6)
print(p_marg)
dev.off()

df2 %>% filter(is_subset == FALSE) %>% wilcox_test(Inclusion_Pct~geneType)
df2 %>% filter(is_subset == FALSE) %>% wilcox_test(Identity_Pct~geneType)
