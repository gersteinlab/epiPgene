rm(list=ls())
library(dplyr)
library(biomaRt)
library(latex2exp)
library(ggplot2)
library(scales)
path <- '/gpfs/gibbs/pi/gerstein/yj329/gnomad/intermediates/upstream'
file_list <- list.files(path, pattern='*.txt', full.names = T)
res <- list()
for (file in file_list) {
  df <- read.delim(file, sep = '\t', header = F)
  tmp <- df %>% dplyr::filter(V8 != '.', V11 != '.') %>% ## have known SNP
    dplyr::mutate(V16 = 1 - as.numeric(V12), V17 = 1 - as.numeric(V13), nfe_maf = pmin(as.numeric(V12), V16), afr_maf = pmin(as.numeric(V13), V17)) %>%
    dplyr::select(-V12, -V13, -V16, -V17) %>%
    dplyr::rename(cadd_phred = V14, sift_max = V15)
  res[[file]] <- tmp
}
res <- do.call(rbind, res, quote = F)
res <- res %>% dplyr::rename(transcript_id = V6, geneType = V7)

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
df <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), filters = "ensembl_transcript_id", values = unique(res$transcript_id), mart = mart)
res <- res %>% merge(df, by.x = "transcript_id", by.y = "ensembl_transcript_id")

res2 <- res %>%
  mutate(maf_bin = case_when(
    nfe_maf < 1e-4 ~ "< 1e-4",
    nfe_maf < 1e-3 ~ "[1e-4, 1e-3)",
    TRUE           ~ "[1e-3, 0.5)"))

maf_prop <- res2 %>% group_by(geneType, maf_bin) %>%
  summarise(n = n()) %>% group_by(geneType) %>% mutate(prop = n / sum(n))

maf_prop <- maf_prop %>%
  mutate(maf_bin = factor(maf_bin, levels = c("< 1e-4", "[1e-4, 1e-3)", "[1e-3, 0.5)")),
         geneType = factor(geneType, levels = c("protein_coding", "lncRNA", "processed_pseudogene", "unprocessed_pseudogene", "unitary_pseudogene")))

## proportion plot
p <- ggplot(maf_prop, aes(x = geneType, y = prop, fill = maf_bin)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = percent(prop, accuracy = 0.01)), position = position_stack(vjust = 0.5), size = 3, color = "black") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_discrete(labels = c(
    protein_coding         = "Protein-coding",
    unprocessed_pseudogene = "Unprocessed",
    processed_pseudogene   = "Processed",
    unitary_pseudogene     = "Unitary",
    lncRNA                 = "lncRNA"
  )) + labs(x = NULL, y = "Proportion of SNVs with different MAF", fill = "MAF bin") +
  scale_fill_manual(values = c('< 1e-4' = '#457b9d', "[1e-4, 1e-3)" = "#a8dadc", "[1e-3, 0.5)" = "#f1faee")) +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color="black", size = 12))

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/plots/upstream_MAF_prop.pdf', p, width = 8, height = 6)


maf_tail <- maf_prop %>% filter(maf_bin == "[1e-3, 0.5)")
p2 <- ggplot(maf_tail, aes(x = geneType, y = prop)) +
  geom_col(width = 0.75, fill = "#f1faee") +
  scale_y_continuous(labels = percent_format(accuracy = 0.01),
                     expand = expansion(mult = c(0, 0.1))) +
  theme_bw(base_size = 14) +
  scale_x_discrete(labels = c(
    protein_coding         = "Protein-coding",
    unprocessed_pseudogene = "Unprocessed",
    processed_pseudogene   = "Processed",
    unitary_pseudogene     = "Unitary",
    lncRNA                 = "lncRNA"
  )) +
  labs(x = NULL, y = "Proportion of SNVs with MAF >= 1e-3") +
  theme_bw(base_size = 14) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_text(color="black", size = 12))
ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig3/plots/upstream_MAF_prop_zoomin.pdf', p2, width = 6, height = 4)

# unprocessed and parent
mapping <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/parent/tss1kb_comparison/tss1kb_inclusion_identity.txt', sep = '\t', header = T)
parentGene_unprocessed <- mapping %>% filter(geneType == 'Unprocessed') %>% pull(ParentGene) %>% unique()
parentGene_processed <- mapping %>% filter(geneType == 'Processed') %>% pull(ParentGene) %>% unique()

cdf_parent_unprocessed <- res2 %>%
  mutate(
    group2 = case_when(
      ensembl_gene_id %in% parentGene_unprocessed ~ "Parent",
      geneType == "unprocessed_pseudogene" ~ "Unprocessed",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(group2)) %>%
  mutate(
    maf_plot = pmax(nfe_maf, 1e-6)  # 避免 log10(0)
  )

ggplot(cdf_parent_unprocessed,
       aes(x = maf_plot, color = group2)) +
  stat_ecdf(size = 1) +
  scale_x_log10(
    breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_bw(base_size = 14) +
  labs(
    x = "Minor allele frequency (log scale)",
    y = "Cumulative fraction of SNVs",
    color = NULL
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black", size = 12),
    legend.position = "top"
  )+ geom_vline(xintercept = c(1e-4, 1e-3),
                linetype = "dashed", color = "grey50")

parent_unprocessed <- res2 %>%
  mutate(group2 = case_when(
    ensembl_gene_id %in% parentGene_unprocessed ~ "Parent",
    geneType == "unprocessed_pseudogene" ~ "Unprocessed",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group2))

maf_prop_parent_unprocessed <- parent_unprocessed  %>%
  group_by(group2, maf_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group2) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    maf_bin = factor(maf_bin, levels = c("< 1e-4", "[1e-4, 1e-3)", "[1e-3, 0.5)")),
    group2  = factor(group2, levels = c("Parent", "Unprocessed"))
  )



tbl_tail <- matrix(
  c(
    6195, 174643 + 10030,   # Parent
    1533, 36998  + 2340    # Unprocessed
  ),
  nrow = 2,
  byrow = TRUE
)

rownames(tbl_tail) <- c("Parent", "Unprocessed")
colnames(tbl_tail) <- c("MAF>=1e-3", "MAF<1e-3")

tbl_tail


# Fisher exact test
fisher.test(tbl_tail)

parent_processed <- res2 %>%
  mutate(group2 = case_when(
    ensembl_gene_id %in% parentGene_processed ~ "Parent",
    geneType == "processed_pseudogene" ~ "Processed",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(group2))


maf_prop_parent_processed <- parent_processed  %>%
  group_by(group2, maf_bin) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group2) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(
    maf_bin = factor(maf_bin, levels = c("< 1e-4", "[1e-4, 1e-3)", "[1e-3, 0.5)")),
    group2  = factor(group2, levels = c("Parent", "Processed")))


tbl_tail <- matrix(
  c(
    19308, 645225 + 33085,   # Parent
    1296, 30801+1957    # Processed
  ),
  nrow = 2,
  byrow = TRUE
)



