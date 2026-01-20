rm(list = ls())
library(dplyr)
transcript_data <- data.frame(
  type = rep(c("Protein-coding", "lncRNA", "Processed", "Unprocessed", "Unitary"), 2),
  count = c(106643, 17624, 9501, 4232, 520, 22156, 1231, 761, 471, 35),
  filter = rep(c("Annotation-based", "Expression-based"), each = 5))

before_data <- transcript_data %>% filter(filter == "Annotation-based")
after_data <- transcript_data %>% filter(filter == "Expression-based")
reduction_data <- before_data %>%
  mutate(
    after = after_data$count,
    reduction_pct = round((count - after) / count * 100, 1),
    label = paste0(reduction_pct, "%"))

transcript_data$type <- factor(transcript_data$type, levels = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed', 'Unitary'))
reduction_data$type <- factor(reduction_data$type, levels = c('Protein-coding', 'lncRNA', 'Processed', 'Unprocessed', 'Unitary'))

library(ggsci)
p <- ggplot(transcript_data, aes(x = filter, y = count, fill = filter)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~type, scales = 'free_y', ncol = 1) +
  geom_text(data = reduction_data, aes(x = "Expression-based", y = count, label = label), vjust = 6, size = 4,inherit.aes = FALSE) +
  labs(x = "", y = "Transcript count (non-redundant)") +
  scale_fill_manual(values = c('#CC99CCFF','#00FFFFFF')) +
  scale_y_continuous(labels = scales::comma) +
  theme_bw(base_size = 14) +
  theme( 
    legend.position = 'None',
    axis.text = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

ggsave('num_tx_filtering_reduced.pdf', p, width = 2, height = 8)
