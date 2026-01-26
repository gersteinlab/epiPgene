cols <- c('Protein-coding' = 'black', 'lncRNA' = 'gray', 'Unprocessed' = '#edae49', 'Processed' = '#d1495b', 'Unitary' = '#00798c')
library(ggplot2)
data <- data.frame(
  Categories = c("Processed", "Protein-coding"),
  Percentages = c(9.9, 1.42))

p1 <- ggplot(data, aes(x = Categories, y = Percentages, fill = Categories)) +
  geom_bar(stat = "identity", width = 0.7) + 
  scale_fill_manual(values = cols) +
  labs(
    title = "YY1(Zf)/Homer(0.696)",
    x = "",
    y = "Percentage (%)"
  ) +
  theme_bw() +
  theme(legend.position = "none",  axis.text = element_text(color="black", size = 12)) +
  ylim(0, 13)

# processed vs unprocessed
data <- data.frame(
  Categories = c("Processed", "Unprocessed"),
  Percentages = c(30.49, 9.35)
)

p2 <- ggplot(data, aes(x = Categories, y = Percentages, fill = Categories)) +
  geom_bar(stat = "identity", width = 0.7) + 
  scale_fill_manual(values = cols) +
  labs(
    title = "YY2/Jaspar(0.784)",
    x = "",
    y = "Percentage (%)"
  ) +
  theme_bw() +
  theme(legend.position = "none",  axis.text = element_text(color="black", size = 12)) +
  ylim(0, 31)

##
data <- data.frame(
  Categories = c("Processed", "Unprocessed"),
  Percentages = c(12.22, 1.47)
)

p3 <- ggplot(data, aes(x = Categories, y = Percentages, fill = Categories)) +
  geom_bar(stat = "identity", width = 0.7) + 
  scale_fill_manual(values = cols) +
  labs(
    title = "Yy1/Jaspar(0.712)",
    x = "",
    y = "Percentage (%)"
  ) +
  theme_bw() +
  theme(legend.position = "none",  axis.text = element_text(color="black", size = 12)) +
  ylim(0, 13)


library(patchwork)
p <- p1 + p3 + p2 + plot_layout(axis_titles= 'collect')
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/plots/homer2_enrichment.pdf', width = 8.5, height = 4)
