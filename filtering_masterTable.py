#/usr/bin/env python
import pandas as pd
import numpy as np
import glob

df = pd.read_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase.txt', sep = '\t')

# extract upstream regions for assessing functional properties
upstream = (
    df.loc[:, ["chr", "start", "end", "gene_id", "score", "strand", "transcript_id", "geneType"]]
    .drop_duplicates()
    .assign(
        new_start=lambda x: np.where(x["strand"] == "-", x["end"] - 1000, x["start"]),
        new_end=lambda x: np.where(x["strand"] == "-", x["end"], x["start"] + 1000)
    )
    .drop(columns=["start", "end"])
    .rename(columns={"new_start": "start", "new_end": "end"})
    .loc[:, ["chr", "start", "end", "gene_id", "score", "strand", "transcript_id", "geneType"]])
upstream.to_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/resources/GRCh38_tss_upstream.bed', sep = '\t', index=False)

# DNA methylation data
methylation_files = glob.glob('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/intersection/DNAme/*.txt')
methylation_data = []
for file in methylation_files:
    tissue = file.split('/')[-1].replace('.txt', '')
    meth_df = pd.read_csv(file, sep='\t')
    meth_df['transcript_id'] = meth_df['transcript_id'].astype(str).str.split('.').str[0]
    meth_df['tissue'] = tissue
    methylation_data.append(meth_df)
methylation_df = pd.concat(methylation_data, ignore_index=True)
tmp = df.merge(methylation_df, on=['transcript_id', 'tissue'], how='left')
tmp.to_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase+DNAme.txt', sep = '\t', index=False)

# transcripts with overlapping promoter regions
blacklist = pd.read_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/gencode_v29/blacklist.txt', header=None)
lt = [item.split('.')[0] for item in blacklist[0] if "PAR_Y" not in item]
# remove transcripts and rename gene type
df2 = df[~df['transcript_id'].isin(lt)]
df2['geneType'] = df2['geneType'].apply(lambda x: 
    'unprocessed_pseudogene' if 'unprocessed_pseudogene' in x else
    'processed_pseudogene' if 'processed_pseudogene' in x else
    'unitary_pseudogene' if 'unitary_pseudogene' in x else
    x)

lncRNA_types = [
    '3prime_overlapping_ncrna', 'antisense', 'bidirectional_promoter_lncrna', 
    'macro_lncRNA', 'non_coding', 'processed_transcript', 
    'sense_intronic', 'sense_overlapping', 'lincRNA'
]
df2['geneType'] = df2['geneType'].replace(lncRNA_types, 'lncRNA')

# annotation-based filtering
df_unique = df2[['transcript_id', 'geneType']].drop_duplicates()
transcript_counts = df_unique.groupby('geneType')['transcript_id'].nunique()
print(transcript_counts)

# expression-based filtering
df2_filtered = df2.loc[df2.groupby(['gene_id', 'tissue'])['avgTPM'].idxmax()].reset_index(drop=True)
df2_filtered = df2_filtered[df2_filtered['avgTPM'] >= 1]

# expression-based filtering -- num of non-redundant transcripts
df_unique = df2_filtered[['transcript_id','geneType']].drop_duplicates()
transcript_counts = df_unique.groupby('geneType')['transcript_id'].nunique()
print(transcript_counts)

# expression-based filtering -- num of non-redundant genes
df_unique = df2_filtered[['gene_id','geneType']].drop_duplicates()
gene_counts = df_unique.groupby('geneType')['gene_id'].nunique()
print(gene_counts)

# (tissue, transcript）-- pair of transcripts and tissues
df_unique = df2_filtered[['tissue', 'transcript_id', 'geneType']].drop_duplicates()
df_unique['tissue_transcript'] = list(zip(df_unique['tissue'], df_unique['transcript_id']))
transcript_counts = df_unique.groupby('geneType')['tissue_transcript'].nunique()
print(transcript_counts)

# RAMPAGE-based filter
RAMPAGE = pd.read_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/intersection/RAMPAGE/GRCh38_tss_1kb_rPeaks.txt', sep='\t', header=None)
RAMPAGE_filtered = RAMPAGE[~RAMPAGE[6].str.contains("PAR_Y", na=False)]
RAMPAGE_filtered[6] = RAMPAGE_filtered[6].str.split('.').str[0]
RAMPAGE_filtered = RAMPAGE_filtered[RAMPAGE_filtered[8] != '.']
tx = set(RAMPAGE_filtered[6])
df2_filtered['RAMPAGE'] = df2_filtered['transcript_id'].isin(tx)

df2_filtered = df2_filtered.sort_values(by=["chr", "start"]).reset_index(drop=True)
df2_filtered.to_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase_filtered.txt', sep = '\t', index=False)

