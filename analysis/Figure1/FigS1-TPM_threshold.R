library(dplyr)
library(ggplot2)

dat <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_ALL_Transcripts_0.1.txt', sep = '\t', header = T)

cutoffs <- c(seq(0.1, 1, 0.1), seq(1.5, 5, 0.5))

th_counts <- map_df(cutoffs, function(t) {
  dat %>%                                 
    filter(avgTPM >= t) %>%  
    distinct(transcript_id, geneType) %>% 
    count(geneType, name = "n_transcripts") %>% 
    mutate(cutoff = t)
})

cols <- 

p <- ggplot(th_counts,
       aes(x = cutoff, y = n_transcripts,
           colour = geneType)) +
  labs(x = 'Average TPM across donors', y = 'Number of transcripts') +
  geom_line(linewidth = 1) +
  geom_point(size = 2) + 
  facet_wrap(~geneType, scales = 'free') +
  theme(panel.border = element_rect(colour = "black", fill = NA), legend.position = 'None')

ggsave('num_transcripts_TPM.pdf', width = 8, height = 8)
