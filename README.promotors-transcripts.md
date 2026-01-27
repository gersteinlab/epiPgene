We have provided a comprehensive table for characterizing promoter activity of GENCODE v29 transcripts across tissues. Please visit [here](http://meetings.gersteinlab.org/2025/02.24/promoters-transcripts.txt.zip) to download.

The dataset contains the following columns, each representing specific information about promoter regions and related biological data:

**1-3.** BED intervals that define the promoter regions (TSS ± 1,000 bp).

**4.** Ensembl gene ID corresponding to the associated transcript.

**5.** Score.

**6.** Strand.

**7.** Ensembl transcript ID to which the promoter was assigned.

**8.** Gene biotype, as annotated in GENCODE v29.

**9.** Tissue name.

**10.** Average transcript-level expression (in TPM) across donors and technical replicates in the tissue.

**11–16.** Length of overlap between the promoter regions and peaks from various histone mark ChIP-seq datasets, including H3K27ac, H3K4me3, H3K27me3, H3K36me3, H3K4me1, and H3K9me3. “NA” indicates missing experimental data, and a value of “0” indicates no overlap between the promoter region and the corresponding peak.

**17–18.** Length of overlap between the promoter regions and peaks from DNase-seq and ATAC-seq data.

**19–22.** Number of CpG sites within the promoters and different levels of methylation.

