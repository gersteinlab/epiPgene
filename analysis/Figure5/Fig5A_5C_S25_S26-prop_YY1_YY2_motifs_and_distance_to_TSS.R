rm(list = ls())
library(ggplot2)
library(dplyr)
library(GenomicRanges)
cols <- c('protein_coding' = 'black', 'lncRNA' = 'gray','unprocessed_pseudogene' = '#edae49', 'processed_pseudogene' = '#d1495b', 'unitary_pseudogene' = '#00798c')

findMotifInPeaks <- function(promoter_file, motif_file, type){
  promoter_df <- read.table(promoter_file, skip = 1, sep = '\t') %>% dplyr::select(V1, V2, V3, V4, V5, V10)
  promoter_gr <- GRanges(seqnames = promoter_df$V2, ranges = IRanges(start = as.numeric(promoter_df$V3)+1, end = as.numeric(promoter_df$V4)), strand = promoter_df$V5, gene_id = promoter_df$V1)
  
  motif_df <- read.table(motif_file, skip = 1, sep = '\t') %>% distinct()
  motif_gr <- GRanges(seqnames = motif_df$V1, ranges = IRanges(start = as.numeric(motif_df$V2)+1, end = as.numeric(motif_df$V3)), strand = motif_df$V6, motif = motif_df$V4, score = motif_df$V5)
  
  overlaps <- findOverlaps(motif_gr, promoter_gr, ignore.strand = TRUE)
  
  res <- data.frame(
    promoter_id = mcols(promoter_gr)$gene_id[subjectHits(overlaps)],
    promoter_chr = as.character(seqnames(promoter_gr))[subjectHits(overlaps)],
    promoter_start = start(promoter_gr)[subjectHits(overlaps)],
    promoter_end = end(promoter_gr)[subjectHits(overlaps)],
    promoter_strand = as.character(strand(promoter_gr))[subjectHits(overlaps)],
    motif_chr = as.character(seqnames(motif_gr))[queryHits(overlaps)],
    motif_start = start(motif_gr)[queryHits(overlaps)],
    motif_end = end(motif_gr)[queryHits(overlaps)],
    motif_strand = as.character(strand(motif_gr))[queryHits(overlaps)],
    motif_seq = mcols(motif_gr)$motif[queryHits(overlaps)],
    motif_score = mcols(motif_gr)$score[queryHits(overlaps)])
  
  res$tss <- round((res$promoter_start + res$promoter_end) / 2)
  res$distance_to_tss <- ifelse(res$promoter_strand == "+", res$motif_start - res$tss, res$tss - res$motif_start)
  
  res <- res %>% mutate(geneType = type)
  return(res)
}
# motif 1
YY1_0.696_processed_pgene <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/processed_pgene/YY1_0.696_processed_pgene_promoters.all.txt',
                                              '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/processed_pgene/YY1_0.696_processed_pgene_promoters.bed', 'processed_pseudogene')

sum(YY1_0.696_processed_pgene$promoter_strand == YY1_0.696_processed_pgene$motif_strand)/nrow(YY1_0.696_processed_pgene)

YY1_0.696_unprocessed_pgene <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/unprocessed_pgene/YY1_0.696_unprocessed_pgene_promoters.all.txt',
                                                '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/unprocessed_pgene/YY1_0.696_unprocessed_pgene_promoters.bed', 'unprocessed_pseudogene')
sum(YY1_0.696_unprocessed_pgene$promoter_strand == YY1_0.696_unprocessed_pgene$motif_strand)/nrow(YY1_0.696_unprocessed_pgene)


YY1_0.696_pcg <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/pcg/YY1_0.696_protein_coding_promoters.all.txt',
                                  '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/pcg/YY1_0.696_protein_coding_promoters.bed', 'protein_coding')

YY1_0.696_lncRNA <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/lncRNA/YY1_0.696_lncRNA_promoters.all.txt',
                                     '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/lncRNA/YY1_0.696_lncRNA_promoters.bed', 'lncRNA')

YY1_0.696 <- rbind(YY1_0.696_pcg, YY1_0.696_lncRNA, YY1_0.696_processed_pgene, YY1_0.696_unprocessed_pgene) %>% filter(!grepl("-", promoter_id))
write.table(YY1_0.696, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/data/YY1_0.696_near_TSS.txt')

p <- rbind(YY1_0.696_processed_pgene, YY1_0.696_pcg) %>% filter(!grepl("-", promoter_id)) %>%
  ggplot(aes(x = distance_to_tss)) + 
  geom_density(aes(color = geneType), size = 0.5, adjust = 1) +
  scale_x_continuous(name = "Distance to TSS (bp)", breaks = seq(-100, 100, by = 10)) + 
  scale_y_continuous(name = "Density") + labs(title = 'YY1_0.696') +
  scale_color_manual(values = cols, name = '',
                     labels = c("processed_pseudogene" = "Processed", "protein_coding" = "Protein-cod.")) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_line(size = 0.5),
    axis.text = element_text(color="black", size = 12))

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/plots/distance_to_tss_YY1_0.696.pdf', p, width = 10, height = 10)

# motif 2
Yy1_0.712_processed_pgene <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/processed_pgene/Yy1_0.712_processed_pgene_promoters.all.txt',
                                              '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/processed_pgene/Yy1_0.712_processed_pgene_promoters.bed', 'processed_pseudogene')
sum(Yy1_0.712_processed_pgene$promoter_strand == Yy1_0.712_processed_pgene$motif_strand)/nrow(Yy1_0.712_processed_pgene)

Yy1_0.712_unprocessed_pgene <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/unprocessed_pgene/Yy1_0.712_unprocessed_pgene_promoters.all.txt',
                                                '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/unprocessed_pgene/Yy1_0.712_unprocessed_pgene_promoters.bed', 'unprocessed_pseudogene')
sum(Yy1_0.712_unprocessed_pgene$promoter_strand == Yy1_0.712_unprocessed_pgene$motif_strand)/nrow(Yy1_0.712_unprocessed_pgene)

Yy1_0.712_pcg <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/pcg/Yy1_0.712_protein_coding_promoters.all.txt',
                                 '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/pcg/Yy1_0.712_protein_coding_promoters.bed', 'protein_coding')
Yy1_0.712_lncRNA <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/lncRNA/Yy1_0.712_lncRNA_promoters.all.txt',
                                  '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/lncRNA/Yy1_0.712_lncRNA_promoters.bed', 'lncRNA')

Yy1_0.712 <- rbind(Yy1_0.712_pcg, Yy1_0.712_lncRNA, Yy1_0.712_processed_pgene, Yy1_0.712_unprocessed_pgene) %>% filter(!grepl("-", promoter_id))
write.table(Yy1_0.712, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/data/Yy1_0.712_near_TSS.txt')

p2 <- rbind(Yy1_0.712_processed_pgene, Yy1_0.712_unprocessed_pgene) %>% filter(!grepl("-", promoter_id)) %>%
  ggplot(aes(x = distance_to_tss)) + 
  geom_density(aes(color = geneType), size = 0.5, adjust = 1) +
  scale_x_continuous(name = "Distance to TSS (bp)", breaks = seq(-100, 100, by = 10)) + 
  scale_y_continuous(name = "Density") + labs(title = 'Yy1_0.712') +
  scale_color_manual(values = cols, name = '',
                     labels = c("processed_pseudogene" = "Processed", "unprocessed_pseudogene" = "Unprocessed")) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_line(size = 0.5),
    axis.text = element_text(color="black", size = 12))
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/plots/distance_to_tss_Yy1_0.712.pdf', p2, width = 10, height = 10)

# motif 3
YY2_0.784_processed_pgene <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/processed_pgene/YY2_0.784_processed_pgene_promoters.all.txt',
                                              '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/processed_pgene/YY2_0.784_processed_pgene_promoters.bed', 'processed_pseudogene')
sum(YY2_0.784_processed_pgene$promoter_strand == YY2_0.784_processed_pgene$motif_strand)/nrow(YY2_0.784_processed_pgene)

YY2_0.784_unprocessed_pgene <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/unprocessed_pgene/YY2_0.784_unprocessed_pgene_promoters.all.txt',
                                                '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/unprocessed_pgene/YY2_0.784_unprocessed_pgene_promoters.bed', 'unprocessed_pseudogene')
sum(YY2_0.784_unprocessed_pgene$promoter_strand == YY2_0.784_unprocessed_pgene$motif_strand)/nrow(YY2_0.784_unprocessed_pgene)

YY2_0.784_pcg <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/pcg/YY2_0.784_protein_coding_promoters.all.txt',
                                  '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/pcg/YY2_0.784_protein_coding_promoters.bed', 'protein_coding')

YY2_0.784_lncRNA <- findMotifInPeaks('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/lncRNA/YY2_0.784_lncRNA_promoters.all.txt',
                                  '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/lncRNA/YY2_0.784_lncRNA_promoters.bed', 'lncRNA')


YY2_0.784 <- rbind(YY2_0.784_pcg, YY2_0.784_lncRNA, YY2_0.784_processed_pgene, YY2_0.784_unprocessed_pgene) %>% filter(!grepl("-", promoter_id))
write.table(YY2_0.784, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/data/YY2_0.784_near_TSS.txt')

p3 <- rbind(YY2_0.784_processed_pgene, YY2_0.784_unprocessed_pgene) %>% filter(!grepl("-", promoter_id)) %>%
  ggplot(aes(x = distance_to_tss)) + 
  geom_density(aes(color = geneType), size = 0.5, adjust = 1) +
  scale_x_continuous(name = "Distance to TSS (bp)", breaks = seq(-100, 100, by = 10)) + 
  scale_y_continuous(name = "Density") + labs(title = 'YY2_0.784') +
  scale_color_manual(values = cols, name = '',
                     labels = c("processed_pseudogene" = "Processed", "unprocessed_pseudogene" = "Unprocessed")) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_line(size = 0.5),
    axis.text = element_text(color="black", size = 12))
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/plots/distance_to_tss_YY2_0.784.pdf', p3, width = 10, height = 10)

# strand bias
df <- rbind(YY1_0.696, Yy1_0.712)
df %>%
  filter(geneType == 'processed_pseudogene') %>%
  summarize(
    n_same = sum(promoter_strand == motif_strand, na.rm = TRUE),
    total  = n(),
    proportion = mean(promoter_strand == motif_strand, na.rm = TRUE)
  )

# calculate p-value
calculate_enrichment_p <- function(df, window = 10) {
  df$in_center <- abs(df$distance_to_tss) < window
  tab <- table(df$geneType, df$in_center)
  fisher <- fisher.test(tab)

  return(list(
    p.value = fisher$p.value,
    odds.ratio = fisher$estimate,
    contingency.table = tab
  ))
}

calculate_enrichment_p(rbind(YY1_0.696_pcg, YY1_0.696_processed_pgene))
calculate_enrichment_p(rbind(Yy1_0.712_unprocessed_pgene, Yy1_0.712_processed_pgene))
calculate_enrichment_p(rbind(YY2_0.784_unprocessed_pgene, YY2_0.784_processed_pgene))


YY1_0.696_processed_pgene$promoter_id
Yy1_0.712_processed_pgene$promoter_id
YY2_0.784_processed_pgene$promoter_id

library(ggVennDiagram)
venn_input <- list(
  YY1_0.696 = YY1_0.696_processed_pgene$promoter_id,
  Yy1_0.712 = Yy1_0.712_processed_pgene$promoter_id,
  YY2_0.784 = YY2_0.784_processed_pgene$promoter_id)

p4 <- ggVennDiagram(venn_input, label_alpha = 0) +
  scale_fill_gradient(low = "#c9def4", high = "#b190ba")
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/plots/processed_pgene_overlap.pdf', p4, width = 8, height = 8)


# calculate the proportion
numerator_processed_pgene <- unique(c(YY1_0.696_processed_pgene$promoter_id, Yy1_0.712_processed_pgene$promoter_id))
denominator_processed_pgene <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/processed_pgene/YY1_0.696_processed_pgene_promoters.all.txt', skip = 1, sep = '\t') %>% pull(V1)

numerator_unprocessed_pgene <- unique(c(YY1_0.696_unprocessed_pgene$promoter_id, Yy1_0.712_unprocessed_pgene$promoter_id))
denominator_unprocessed_pgene <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/unprocessed_pgene/Yy1_0.712_unprocessed_pgene_promoters.all.txt', skip = 1, sep = '\t') %>% pull(V1)

numerator_pcg <- unique(c(YY1_0.696_pcg$promoter_id, Yy1_0.712_pcg$promoter_id))
denominator_pcg <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/pcg/YY1_0.696_protein_coding_promoters.all.txt', skip = 1, sep = '\t') %>% pull(V1)

numerator_lncRNA <- unique(c(YY1_0.696_lncRNA$promoter_id, Yy1_0.712_lncRNA$promoter_id))
denominator_lncRNA <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/lncRNA/YY1_0.696_lncRNA_promoters.all.txt', skip = 1, sep = '\t') %>% pull(V1)

df <- tibble(
  category = c("processed_pseudogene", "unprocessed_pseudogene", "protein_coding", 'lncRNA'),
  numerator = c(length(numerator_processed_pgene), length(numerator_unprocessed_pgene), length(numerator_pcg), length(numerator_lncRNA)),
  denominator = c(length(denominator_processed_pgene), length(denominator_unprocessed_pgene), length(denominator_pcg), length(denominator_lncRNA))
) %>%
  mutate(proportion = numerator / denominator)

df$category <- factor(df$category, levels =c('processed_pseudogene', 'lncRNA', 'protein_coding', 'unprocessed_pseudogene'))
library(scales)
p5 <- ggplot(df, aes(x = category, y = proportion, fill = category)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(labels = c("processed_pseudogene" = "Processed", "unprocessed_pseudogene" = "Unprocessed", "protein_coding" = 'Protein-cod.', "lncRNA" = "lncRNA")) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal() +
  labs(x = '', y = 'Proportion', title = 'YY1-like motifs') +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_line(size = 0.5),
    axis.text = element_text(color="black", size = 12),
    legend.position = 'None')
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/plots/propotion_of_promoters_YY1_like_motifs.pdf', p5, width = 6, height = 6)

# YY2-like
numerator_processed_pgene <- unique(c(YY2_0.784_processed_pgene$promoter_id))
denominator_processed_pgene <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/processed_pgene/YY2_0.784_processed_pgene_promoters.all.txt', skip = 1, sep = '\t') %>% pull(V1)

numerator_unprocessed_pgene <- unique(c(YY2_0.784_unprocessed_pgene$promoter_id))
denominator_unprocessed_pgene <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/unprocessed_pgene/YY2_0.784_unprocessed_pgene_promoters.all.txt', skip = 1, sep = '\t') %>% pull(V1)

numerator_pcg <- unique(c(YY2_0.784_pcg$promoter_id))
denominator_pcg <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/pcg/YY2_0.784_protein_coding_promoters.all.txt', skip = 1, sep = '\t') %>% pull(V1)

numerator_lncRNA <- unique(c(YY2_0.784_lncRNA$promoter_id))
denominator_lncRNA <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/motif_analysis/lncRNA/YY2_0.784_lncRNA_promoters.all.txt', skip = 1, sep = '\t') %>% pull(V1)

df <- tibble(
  category = c("processed_pseudogene", "unprocessed_pseudogene", "protein_coding", 'lncRNA'),
  numerator = c(length(numerator_processed_pgene), length(numerator_unprocessed_pgene), length(numerator_pcg), length(numerator_lncRNA)),
  denominator = c(length(denominator_processed_pgene), length(denominator_unprocessed_pgene), length(denominator_pcg), length(denominator_lncRNA))
) %>%
  mutate(proportion = numerator / denominator)

df$category <- factor(df$category, levels =c('processed_pseudogene', 'protein_coding', 'lncRNA', 'unprocessed_pseudogene'))
library(scales)

p6 <- ggplot(df, aes(x = category, y = proportion, fill = category)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  geom_text(aes(label = scales::percent(proportion, accuracy = 0.1)), vjust = -0.5, size = 4) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(labels = c("processed_pseudogene" = "Processed", "unprocessed_pseudogene" = "Unprocessed", "protein_coding" = 'Protein-cod.', "lncRNA" = "lncRNA")) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal() +
  labs(x = '', y = 'Proportion', title = 'YY2-like motif') +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),  
    panel.grid.major.y = element_line(size = 0.5, linetype = 'dotted', colour = "gray"),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_line(size = 0.5),
    axis.text = element_text(color="black", size = 12),
    legend.position = 'None')
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/plots/propotion_of_promoters_YY2_like_motifs.pdf', p6, width = 10, height = 10)

## using processed pgene as an example, to do the fisher test
g <- 'processed_pseudogene'
hic_df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig5/tables/hic_coords_transcripts.txt')
all_df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', header = T)
tx_name_with_hic <- hic_df %>% filter(V13 == g) %>% dplyr::select(V1, V2, V13) %>% distinct()
tx_name_without_hic <- all_df %>% filter(geneType == g) %>% dplyr::select(gene_id, transcript_id, geneType) %>% distinct()
tx_with_motif <- rbind(YY1_0.696_processed_pgene, Yy1_0.712_processed_pgene) %>% dplyr::select(promoter_id) %>% distinct()

with_hic_ids    <- tx_name_with_hic$V1
without_hic_ids <- setdiff(tx_name_without_hic$gene_id, with_hic_ids)
motif_ids       <- tx_with_motif$promoter_id


n_with_hic_and_motif     <- sum(with_hic_ids %in% motif_ids)
n_without_hic_and_motif  <- sum(without_hic_ids %in% motif_ids)

n_with_hic_no_motif      <- length(with_hic_ids) - n_with_hic_and_motif
n_without_hic_no_motif   <- length(without_hic_ids) - n_without_hic_and_motif

mat <-matrix(
  c(n_with_hic_and_motif, n_with_hic_no_motif,
    n_without_hic_and_motif, n_without_hic_no_motif),
  nrow = 2, byrow = TRUE
)
colnames(mat) <- c("motif", "no_motif")
rownames(mat) <- c("with_hic", "without_hic")

fisher.test(mat)
mat
