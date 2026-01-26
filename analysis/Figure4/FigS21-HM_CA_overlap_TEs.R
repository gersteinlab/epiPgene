rm(list = ls())
# path
overlap_dir <- "histoneMarks/TE"
bed_dir <- "/gpfs/gibbs/pi/gerstein/yj329/epiPgene/chromatin/assay_tissue.bed"

files <- list.files(overlap_dir, pattern = "*.txt", full.names = TRUE)
results <- data.frame()

for (file in files) {
  filename <- basename(file)
  parts <- strsplit(sub(".txt", "", filename), "\\.")[[1]]
  mark <- parts[1]
  tissue <- parts[2]
  # num of overlap
  df_overlap <- read.table(file, sep = '\t', header = FALSE)
  df_overlap$peak_id <- paste(df_overlap$V1, df_overlap$V2, df_overlap$V3, sep = "_")
  overlap_peak_count <- length(unique(df_overlap$peak_id))
  
  # all peaks
  bed_file <- paste0(bed_dir, "/", gsub('.txt', '.bed', filename))
  df_total <- read.table(bed_file, sep = '\t', header = FALSE)
  df_total$peak_id <- paste(df_total$V1, df_total$V2, df_total$V3, sep = "_")
  total_peak_count <- length(unique(df_total$peak_id))
  
  results <- rbind(results, data.frame(
    assay = mark,
    tissue = tissue,
    overlap_peak_count = overlap_peak_count,
    total_peak_count = total_peak_count,
    overlap_ratio = round(overlap_peak_count/total_peak_count, 4)
  ))
}

library(ggplot2)
library(ggbeeswarm)
palette <- c("ATAC" = "#AF4D85",
             "DNase" = "#C71585",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#D199B9",
             "H3K27me3" = "#1D2976",
             "H3K9me3" = "#A7ADD4",
             "H3K36me3" = "#7FBC41",
             "H4K20me1" = "#4C7027",
             "H3K4me1" = "#E5AB00")


p <- ggplot(results, aes(x = assay, y = overlap_ratio)) +
  geom_violin(aes(fill = assay), alpha = 0.8, width = 1, scale = "width", trim = FALSE) +
  geom_beeswarm(cex = 1, size = 1, priority = "descending")+
  labs(x = '', y = 'Proportion of peaks overlapping TEs (based on hit count)') +
  scale_fill_manual(values = palette) + ylim(c(0,1)) +
  theme_bw() +
  theme(axis.text = element_text(color = 'black'), legend.position = 'None')

ggsave('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/analysis/fig4/plots/HM_peaks_overlap_TEs.pdf', p, width = 6, height = 4)
