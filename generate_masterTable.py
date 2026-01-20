#/usr/bin/env python
import pandas as pd
import numpy as np
import glob, os
from tqdm import tqdm

# expression
expression_files = glob.glob("/gpfs/gibbs/pi/gerstein/yj329/epiPgene/total_RNA-seq/tissue/*.tsv")
expr_data = pd.DataFrame()
for file in expression_files:
    tissue_name = os.path.basename(file).replace(".tsv", "")
    tmp = pd.read_csv(file, sep="\t")
    tmp["tissue"] = tissue_name
    tmp = tmp[["transcript_id", "gene_id", "tissue","avgTPM"]]
    expr_data = pd.concat([expr_data, tmp], axis=0, ignore_index=True)

# histone
histone_files = glob.glob("/gpfs/gibbs/pi/gerstein/yj329/epiPgene/intersection/histoneMarks/*.txt")
histone_dict = {}
for file in tqdm(histone_files):
    histone_name, tissue_name = os.path.basename(file).split('.')[0], os.path.basename(file).split('.')[1]
    print(f'{tissue_name}\t{histone_name}')
    
    tmp = pd.read_csv(file, sep='\t', header=None)
    tmp = tmp[~tmp[3].str.contains("PAR_Y", na=False)]
    tmp[3] = tmp[3].str.split('.').str[0]
    tmp[6] = tmp[6].str.split('.').str[0]
    tmp = tmp.groupby(list(range(8)), as_index=False)[12].sum()
    tmp.columns = ['chr', 'start', 'end', 'gene_id', 'score', 'strand', 'transcript_id', 'geneType'] + [histone_name]

    if tissue_name not in histone_dict:
        histone_dict[tissue_name] = tmp
    else:
        histone_dict[tissue_name] = pd.merge(histone_dict[tissue_name], tmp, on=['chr', 'start', 'end', 'gene_id', 'score', 'strand', 'transcript_id', 'geneType'], how='outer')

# combine histone dict
df_list = []
for tissue, df in histone_dict.items():
    df["tissue"] = tissue
    df_list.append(df)
master_table = pd.concat(df_list, ignore_index=True)

# merge two tables
common_tissues = common_tissues = set(master_table["tissue"]).intersection(set(expr_data["tissue"]))
master_table = master_table[master_table["tissue"].isin(common_tissues)] # filter on tissues

res = pd.merge(master_table, expr_data, on=["transcript_id", 'gene_id', 'tissue'], how="left")
res = res[['chr', 'start', 'end', 'gene_id', 'score', 'strand', 'transcript_id', 'geneType', 'tissue', 'avgTPM', 'H3K27ac', 'H3K4me3', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K9me3', 'ATAC', 'DNase']]
# only consider protein-coding, lncRNA and pseudogene
GENE_TYPE = ['processed_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unitary_pseudogene', 'transcribed_unprocessed_pseudogene',
             'translated_processed_pseudogene', 'translated_unprocessed_pseudogene', 'unitary_pseudogene', 'unprocessed_pseudogene', 
             'protein_coding', 
             '3prime_overlapping_ncrna', 'antisense', 'bidirectional_promoter_lncrna', 'macro_lncRNA', 'non_coding', 'processed_transcript', 'sense_intronic', 'sense_overlapping','lincRNA']
res = res[res["geneType"].isin(GENE_TYPE)]

res = res.sort_values(by=["chr", "start"]).reset_index(drop=True)
res.to_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_RNAseq+histone+ATAC+DNase.txt', sep = '\t', index=False)

# binarize peaks
binary_columns = ['H3K27ac', 'H3K4me3', 'H3K27me3', 'H3K36me3', 'H3K4me1', 'H3K9me3', 'ATAC', 'DNase']
res[binary_columns] = res[binary_columns].applymap(lambda x: 0 if x == 0 else (1 if pd.notna(x) else np.nan)).astype("Int64")
res.to_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/freeze_mastertable_binary_RNAseq+histone+ATAC+DNase.txt', sep = '\t', index=False)
