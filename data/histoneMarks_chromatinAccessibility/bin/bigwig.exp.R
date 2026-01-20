# 1. set wd
setwd("/gpfs/gibbs/pi/gerstein/yj329/epi/chromatin/preprocessing_Jul")


# 2. read list of bigbed files
m <- read.delim("bigwig.files.tsv", h=F)

# 3. keep only chipseq / open chr assays
m <- m[m$V3 == "Histone_ChIP-seq" | m$V5 == "openChr", ]
m$V5 <- gsub("-human", "", m$V5)

# 4. change entex tissue names to gtex tissue names
entex.gtex <- read.delim("GTEX.EN-TEx.tissues.correspondence.tsv")
z <- entex.gtex$GTEx
names(z) <- entex.gtex$EN.TEx
m$V7 <- z[m$V4]

# 5. save output
write.table(m[complete.cases(m), c(1,3,5,7)], "bigwig.exp.tsv", row.names = F, col.names = F, sep="\t", quote = F)
