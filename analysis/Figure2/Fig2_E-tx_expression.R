rm(list = ls())
library(ggplot2)
library(dplyr)
library(stringr)
library(rstatix)
library(ggpubr)
cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')
df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)

data_to_plot <- df %>% filter(RAMPAGE == "True", avgTPM >= 1) %>%
  dplyr::select(geneType, gene_id, transcript_id, tissue, avgTPM) %>% distinct() %>% mutate(expr = log1p(avgTPM))

pval_tbl <- data_to_plot %>%
  wilcox_test(expr~geneType) %>% 
  adjust_pvalue(method = 'BH') %>%
  add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'n.s.')) 

eff_tbl <- data_to_plot %>%
  wilcox_effsize(expr~geneType) 

outliers <- data_to_plot %>%
  group_by(geneType) %>%
  mutate(
    Q1 = quantile(expr, 0.25),
    Q3 = quantile(expr, 0.75),
    IQR = Q3 - Q1,
    is_outlier = expr < (Q1 - 1.5 * IQR) | expr > (Q3 + 1.5 * IQR)) %>% filter(is_outlier)

p <- ggplot(data_to_plot, aes(x = geneType, y = expr))+ geom_violin(aes(fill = geneType), size = 0.3, alpha = 0.6, trim = F) +
  geom_boxplot(width = 0.1, outlier.shape = NA, size = 0.3) +
  scale_x_discrete(limits = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene'), 
                   labels = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed', 'Unitary')) +
  scale_y_continuous(breaks = seq(1, 12,2)) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  stat_compare_means(comparisons = list(c("processed_pseudogene", "protein_coding"), c("processed_pseudogene", "lncRNA"), c("processed_pseudogene", "unprocessed_pseudogene"), c("processed_pseudogene", "unitary_pseudogene")),
                     method = 'wilcox.test', label = 'p.signif', tip.length = 0.01, step.increase = 0.1,
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'n.s.'))) +
  labs(x = '', y = 'Log-transformed transcript expression', fill = 'Group') + theme_classic()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'None',
    axis.text = element_text(color="black", size = 12))

write.table(pval_tbl, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig2/expr_pvalue_Expr.txt', row.names = F, quote = F)
write.table(eff_tbl, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig2/effsize_Expr.txt', row.names = F, quote = F)
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig2/transcript_expr_filtered_by_RAMPAGE.pdf', p, width = 6, height = 4)