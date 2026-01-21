rm(list=ls())
library(dplyr)
library(readr)
library(purrr)
library(ggplot2)
library(ggpubr)

df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
pgene <- df %>% filter(geneType == 'processed_pseudogene' | geneType == 'unprocessed_pseudogene') %>% select(transcript_id, geneType) %>% distinct()
mapping <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/parent/pgene_parent_mapping_Human103.txt', header = T) %>%
  mutate(PgeneEnsembl103ID = sub("\\..*", "", PgeneEnsembl103ID)) %>%
  merge(pgene, by.y = 'transcript_id', by.x = 'PgeneEnsembl103ID') 

expr_dir <- "/gpfs/gibbs/pi/gerstein/yj329/epiPgene/total_RNA-seq/tissue"
expr_files <- list.files(expr_dir, pattern = "\\.tsv$", full.names = TRUE)

expr_list <- lapply(expr_files, function(file) {
  df <- read_tsv(file, col_types = cols())
  tissue <- tools::file_path_sans_ext(basename(file)) 
  
  df <- df %>%
    select(gene_id, transcript_id, avgTPM) %>%
    rename(!!tissue := avgTPM)  
  
  return(df)
})

expr_matrix <- reduce(expr_list, full_join, by = c("gene_id", "transcript_id"))
cor_results <- data.frame(
  pgene_id = character(),
  parent_gene_id = character(),
  correlation = numeric(),
  geneType = character(),
  stringsAsFactors = FALSE)

for (i in 1:nrow(mapping)) {
  pgene_tx <- mapping$PgeneEnsembl103ID[i]
  parent_gene <- mapping$ParentGene[i]
  geneType <- mapping$geneType[i]
  
  pgene_expr <- expr_matrix %>%
    filter(transcript_id == pgene_tx) %>%
    select(-gene_id, -transcript_id)
  
  parent_exprs <- expr_matrix %>%
    filter(gene_id == parent_gene) %>%
    select(-gene_id, -transcript_id)
  
  if (nrow(pgene_expr) == 1 && nrow(parent_exprs) > 0){
    parent_avg_expr <- colMeans(parent_exprs, na.rm = TRUE)
    if (sd(as.numeric(pgene_expr), na.rm = TRUE) == 0 || sd(parent_avg_expr, na.rm = TRUE) == 0) {
      cat("STD = 0 - Skipped\n")
      cat("Pseudogene transcript:", pgene_tx, "\n")
      cat("Parent gene:", parent_gene, "\n")
      cat("pgene_expr:", as.numeric(pgene_expr), "\n")
      cat("parent_avg_expr:", parent_avg_expr, "\n\n")
      next
    }
    
    cor_val <- cor(as.numeric(pgene_expr), parent_avg_expr, use = "pairwise.complete.obs", method = 'spearman')
    
    cor_results <- rbind(cor_results, data.frame(
      pgene_id = pgene_tx,
      parent_gene_id = parent_gene,
      correlation = cor_val,
      geneType = geneType
    ))
  }
}
cols <- c('unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b')

p <- ggplot(cor_results, aes(x = geneType, y = correlation, fill = geneType)) + 
  geom_boxplot(width = 0.2) + labs(x='', y = 'Coexpression correlation') +
  scale_fill_manual(values = cols) + ylim(c(-1,1)) +
  stat_compare_means(method = "wilcox.test") +
  scale_x_discrete(labels = c("processed_pseudogene" = "Processed", "unprocessed_pseudogene" = "Unprocessed")) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1), 
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(color="black", size = 12),
    legend.position = 'None')

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig2/pgene_parent_correlation.pdf', p, width = 4, height = 4)
wilcox.test(correlation ~ geneType, data = cor_results)