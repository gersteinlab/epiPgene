library(ggplot2)
library(ggridges)
library(dplyr)
library(patchwork)

df  <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
figures <- list()
cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray', 'unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')
for (TISSUE in c('Colon_Transverse', 'Muscle_Skeletal')) {
  AB_compartment <- read.table(paste0('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/in_situ_HiC/bwAverageOverBed/', TISSUE, '.bed'))
  colnames(AB_compartment) <- c('chr', 'start', 'end', 'transcript_id', 'signal')
  
  meansd <- function(x, ...) {
    mean <- mean(x)
    sd <- sd(x)
    c(mean - sd, mean, mean + sd)
  }
  
  p <- df %>% filter(tissue == TISSUE, RAMPAGE == "True") %>%
    merge(AB_compartment, by = c('chr', 'start', 'end', 'transcript_id')) %>%
    group_by(geneType) %>%
    filter(n() >= 5) %>%
    ggplot(aes(x = signal, y = geneType, fill = geneType, height = after_stat(density))) + 
    stat_density_ridges(scale = 0.9, quantile_lines = T, quantile_fun = meansd, alpha = 0.6) +
    scale_y_discrete(limits = c('protein_coding', 'lncRNA', 'processed_pseudogene', 'unprocessed_pseudogene'),
                     labels = c('protein_coding' = 'Protein-coding', 'lncRNA' = 'lncRNA', 'processed_pseudogene' = 'Processed', 'unprocessed_pseudogene'='Unprocessed'),
                     expand = c(0.01, 0)) +
    scale_x_continuous(expand = c(0.01, 0)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#023047", size = 0.5) +
    labs(x = 'A/B genome compartments', y = '', title = TISSUE) + 
    scale_fill_manual(values = cols) +
    theme_ridges() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = 'None',
      axis.text = element_text(color="black", size = 12))
 
  figures[[TISSUE]] <- p 
}

(figures[["Colon_Transverse"]] / figures[["Muscle_Skeletal"]]) + plot_layout(axes = "collect", axis_titles = "collect")
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig2/AB_compartments.pdf', width = 5, height = 6)

## p-value
TISSUE <- 'Muscle_Skeletal'
AB_compartment <- read.table(paste0('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/in_situ_HiC/bwAverageOverBed/', TISSUE, '.bed'))
colnames(AB_compartment) <- c('chr', 'start', 'end', 'transcript_id', 'signal')

df %>% filter(tissue == TISSUE, RAMPAGE == "True") %>%
  merge(AB_compartment, by = c('chr', 'start', 'end', 'transcript_id')) %>%
  group_by(geneType) %>%
  filter(n() >= 5) %>% ungroup() %>%
  rstatix::wilcox_test(signal~geneType, ref.group = 'processed_pseudogene')
  