rm(list = ls())
setwd('/gpfs/gibbs/pi/gerstein/yj329/epiPgene')
## YY1_0.696
df1 <- read.table('parent/motif_enrichment/YY1_0.696_parent_processed_pgenes.all.txt', sep = '\t', skip = 1, header = F)
df2 <- read.table('parent/motif_enrichment/YY1_0.696_parent_processed_pgenes.bed', skip = 1, header = F)
##
df3 <- read.table('parent/motif_enrichment/Yy1_0.712_parent_processed_pgenes.all.txt', sep = '\t', skip = 1, header = F)
df4 <- read.table('parent/motif_enrichment/Yy1_0.712_parent_processed_pgenes.bed', sep = '\t', skip = 1, header = F)

prop_1 <- nrow(df2) / nrow(df1)
prop_2 <- nrow(df4) / nrow(df3)

df_plot <- data.frame(
  motif = c("Human-derived", "Mouse-derived"),
  proportion = c(prop_1, prop_2)
)

df_plot$proportion_pct <- df_plot$proportion * 100

library(ggsci)
p <- ggplot(df_plot, aes(x = motif, y = proportion_pct, fill = motif)) +
  geom_bar(stat = "identity", width = 0.6) +  
  geom_text(aes(label = paste0(round(proportion_pct, 1), "%")), 
            vjust = -0.5, size = 5) +
  scale_fill_igv() + ylim(c(0,3)) +
  theme_pubclean(base_size = 14) + labs(x = "", y = "Proportion (%)") +
  theme(legend.position = "none") +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(color="black", size = 12),
    legend.position = 'None')

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/plots/YY1_motif_parent_processedPgene.pdf', p, width = 4, height = 5)
