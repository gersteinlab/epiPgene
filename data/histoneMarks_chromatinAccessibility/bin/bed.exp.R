# 1. set wd
setwd("/gpfs/gibbs/pi/gerstein/yj329/epi/chromatin/preprocessing")

# 2. read list of bigbed files
m <- read.delim("bed.files.tsv", h=F)

# 3. keep only chipseq / open chr assays
m <- m[m$V3 == "Histone_ChIP-seq" | m$V5 == "openChr", ]
m$V5 <- gsub("-human", "", m$V5)
m$V3 <- gsub("-seq", "", m$V3)
m$V5 <- ifelse(m$V5 == "openChr", m$V3, m$V5)

# 4. change entex tissue names to gtex tissue names
entex.gtex <- read.delim("GTEX.EN-TEx.tissues.correspondence.tsv")
z <- entex.gtex$GTEx
names(z) <- entex.gtex$EN.TEx
m$V7 <- z[m$V4]

# 5. save output
write.table(m[complete.cases(m), c(1,5,7)], "bed.exp.tsv", row.names = F, col.names = F, sep="\t", quote = F)