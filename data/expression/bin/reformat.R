# 1. set wd
setwd("/gpfs/gibbs/pi/gerstein/yj329/epi/total_RNA-seq/preprocessing")

# 2. read list of bigbed files
m <- read.delim("transcript_quantification.files.tsv", h=F)

# 4. change entex tissue names to gtex tissue names
entex.gtex <- read.delim("GTEX.EN-TEx.tissues.correspondence.tsv")
z <- entex.gtex$GTEx
names(z) <- entex.gtex$EN.TEx
m$V6 <- z[m$V4]

# 5. save output
write.table(m[complete.cases(m), c(1,3,6)], "transcript_quantification.reformatted.tsv", row.names = F, col.names = F, sep="\t", quote = F)
