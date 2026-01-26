rm(list = ls())
df <- read.table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', head = TRUE)

pcg <- df %>% filter(geneType == 'protein_coding') %>% select(chr, start, end, gene_id, score, strand) %>% distinct()
processed <- df %>% filter(geneType == 'processed_pseudogene') %>% select(chr, start, end, gene_id, score, strand) %>% distinct()
unprocessed <- df %>% filter(geneType == 'unprocessed_pseudogene') %>% select(chr, start, end, gene_id, score, strand) %>% distinct()
lncRNA <- df %>% filter(geneType == 'lncRNA') %>% select(chr, start, end, gene_id, score, strand) %>% distinct()


write.table(pcg, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/homer2/protein_coding_promoter_regions.bed', quote = F, col.names = F, row.names = F, sep = '\t')
write.table(processed, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/homer2/processed_pgene_promoter_regions.bed', quote = F, col.names = F, row.names = F, sep = '\t')
write.table(unprocessed, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/homer2/unprocessed_pgene_promoter_regions.bed', quote = F, col.names = F, row.names = F, sep = '\t')
write.table(lncRNA, '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/homer2/lncRNA_promoter_regions.bed', quote = F, col.names = F, row.names = F, sep = '\t')
