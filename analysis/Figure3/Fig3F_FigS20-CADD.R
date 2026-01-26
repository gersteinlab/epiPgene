rm(list = ls())
library(ggplot2)
library(rstatix)
library(dplyr)
library(ggpubr)
library(latex2exp)
cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')
calc_MAF_and_CADD <- function(path, name, frac = 0.005){
  file_list <- list.files(path , pattern = '*.txt', full.names = T)
  res <- list()
  for (file in file_list) {
    df <- read.delim(file, sep = '\t', header = F)
    tmp <- df %>% dplyr::filter(V8 != '.', V11 != '.') %>% ## have known SNP
      dplyr::mutate(V16 = 1 - as.numeric(V12), V17 = 1 - as.numeric(V13), nfe_maf = pmin(as.numeric(V12), V16), afr_maf = pmin(as.numeric(V13), V17)) %>%
      select(-V12, -V13, -V16, -V17) %>%
      dplyr::rename(cadd_phred = V14, sift_max = V15)
    res[[file]] <- tmp
  }
  res <- do.call(rbind, res, quote = F)
  res <- res %>% dplyr::rename(transcript_id = V6, geneType = V7)
  res$geneType <- factor(res$geneType, levels = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene'))
 
#  outliers <- res %>% filter(nfe_maf > 1e-5) %>% mutate(nfe_maf_log = log10(nfe_maf)) %>%
#    group_by(geneType) %>%
#    mutate(
#      Q1 = quantile(nfe_maf_log, 0.25),
#      Q3 = quantile(nfe_maf_log, 0.75),
#      IQR = Q3 - Q1,
#      is_outlier = nfe_maf_log < (Q1 - 1.5 * IQR) | nfe_maf_log > (Q3 + 1.5 * IQR)) %>% 
#    filter(is_outlier) %>%
#    group_by(geneType) %>% sample_frac(frac)
  
#  MAF_boxplot <- res %>% filter(nfe_maf > 1e-5) %>% mutate(nfe_maf_log = log10(nfe_maf)) %>%
#    ggplot(aes(x = geneType, y = -nfe_maf_log)) +
#    geom_boxplot(linetype="dashed", aes(color = geneType), outlier.shape = NA) +
#    stat_boxplot(aes(ymin=..lower..,ymax=..upper.., fill= geneType),alpha = 0.6, notch = F, outlier.shape = NA) +
#    geom_point(data = outliers, aes(color = geneType), position = position_dodge(width = 0.75), alpha = 0.6) +
#    stat_boxplot(geom = "errorbar",aes(ymax=..ymin.., color = geneType), width=0.4) +
#    scale_fill_manual(values = cols) +
#    scale_color_manual(values = cols) +
#    scale_x_discrete(labels = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed', 'Unitary'))+
#    theme(legend.position = 'None') +
#    stat_compare_means(comparisons = list(c("processed_pseudogene", "protein_coding"), c("processed_pseudogene", "lncRNA"), c("processed_pseudogene", "unprocessed_pseudogene"), c("processed_pseudogene", "unitary_pseudogene") ),
#                       method = 'wilcox.test', label = 'p.signif', tip.length = 0.01, step.increase = 0.05,
#                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'n.s.'))) +
#    labs(x = '', y = TeX("$-\\log_{10}(MAF)$")) +
#    theme_bw() +   
#    theme(
#      panel.grid.major = element_blank(),
#      panel.grid.minor = element_blank(),
#      legend.position = 'None',
#      axis.text = element_text(color="black", size = 12))
  
#  ggsave(paste0('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/plots/', name, '_MAF.pdf'), MAF_boxplot, width = 6, height = 6)
#  tbl <- res %>% filter(nfe_maf > 1e-5) %>% wilcox_test(nfe_maf~geneType, ref.group = 'processed_pseudogene', p.adjust.method = 'fdr')
#  write.table(tbl, paste0('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/tables/', name, '_MAF_wilcoxon_test.txt'), sep = '\t', quote = F)
              
  ## CADD score
  outliers <- res %>% filter(nfe_maf > 1e-5) %>% 
    group_by(geneType) %>%
    mutate(
      Q1 = quantile(as.numeric(cadd_phred), 0.25),
      Q3 = quantile(as.numeric(cadd_phred), 0.75),
      IQR = Q3 - Q1,
      is_outlier = as.numeric(cadd_phred)< (Q1 - 1.5 * IQR) | as.numeric(cadd_phred) > (Q3 + 1.5 * IQR)) %>% 
    filter(is_outlier)           
  
  CADD_boxplot <- res %>% filter(nfe_maf > 1e-5) %>%
    ggplot(aes(x = geneType, y = as.numeric(cadd_phred), fill = geneType, color = geneType)) +
    geom_half_violin(position=position_nudge(x=0.12), side='r', alpha = 0.6, trim = T) +
    geom_boxplot(width=.15, position=position_nudge(x=0), alpha = 0.6, outliers = F) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    scale_x_discrete(labels = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed', 'Unitary'))+
    labs(x = '', y = 'Phred-scaled CADD score') + scale_y_continuous(trans = "sqrt", breaks = scales::trans_breaks("sqrt", function(x) x^2, n = 5)) +
    stat_compare_means(comparisons = list(c("processed_pseudogene", "protein_coding"), c("processed_pseudogene", "lncRNA"), c("processed_pseudogene", "unprocessed_pseudogene"), c("processed_pseudogene", "unitary_pseudogene") ),
                       method = 'wilcox.test', label = 'p.signif', tip.length = 0.01, step.increase = 0.05,
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'n.s.'))) +
    theme_bw() +   
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = 'None',
      axis.text = element_text(color="black", size = 12)) 
  
  ggsave(paste0('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/plots/', name, '_CADD.pdf'), CADD_boxplot, width = 6, height = 6)
 # res$cadd_phred <- as.numeric(res$cadd_phred)
#  tbl <- res %>% filter(nfe_maf > 1e-5) %>% wilcox_test(cadd_phred~geneType, ref.group = 'processed_pseudogene', p.adjust.method = 'fdr')
 # write.table(tbl, paste0('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/tables/', name, '_CADD_wilcoxon_test.txt'), sep = '\t', quote = F)
}

calc_MAF_and_CADD('/gpfs/gibbs/pi/gerstein/yj329/gnomad/intermediates/upstream', 'upstream')
calc_MAF_and_CADD('/gpfs/gibbs/pi/gerstein/yj329/gnomad/intermediates/downstream', 'downstream')
#calc_MAF_and_CADD('/gpfs/gibbs/pi/gerstein/yj329/gnomad/intermediates/exon', 'exon')