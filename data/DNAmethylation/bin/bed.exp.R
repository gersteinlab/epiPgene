# 1. set wd
setwd("/gpfs/gibbs/pi/gerstein/yj329/epi/DNAme/preprocessing")

# 2. read list of bigbed files
m <- read.delim("bed.files.tsv", h=F)

# 3. keep only DNAme array assays
m <- m[m$V3 == "DNAme_array", ]

# 4. change entex tissue names to gtex tissue names
entex.gtex <- read.delim("GTEX.EN-TEx.tissues.correspondence.tsv")
z <- entex.gtex$GTEx
names(z) <- entex.gtex$EN.TEx
m$V6 <- z[m$V4]

# 5. save output
write.table(m[complete.cases(m), c(1,3,6)], "bed.exp.tsv", row.names = F, col.names = F, sep="\t", quote = F)
