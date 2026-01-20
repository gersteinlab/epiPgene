#!/usr/bin/env python

import pandas as pd
import numpy as np
from tqdm import tqdm
import os

metadata = pd.read_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/total_RNA-seq/preprocessing/transcript_quantification.reformatted.tsv', sep = '\t', header=None)
metadata.columns = ['id', 'assay', 'tissue']

gene_type_v29 = pd.read_table('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/gencode_v29/geneType.v29.txt', header = None)
gene_type_v29.columns = ['gene_id', 'transcript_id', 'gene_type_v29']
gene_type_v29['gene_id'] = gene_type_v29['gene_id'].str.replace(r'\.\d+', '', regex=True) # stable gene id
gene_type_v29['transcript_id'] = gene_type_v29['transcript_id'].str.replace(r'\.\d+', '', regex=True) # stable tx id

# read all transcript-level quantification files (in tsv)
root = '/gpfs/gibbs/pi/gerstein/yj329/epiPgene/total_RNA-seq/preprocessing/transcript_quantification/'
df = pd.DataFrame(columns=['transcript_id', 'gene_id', 'length']) # using rsem for quantification

for file in tqdm(os.listdir(root)):
    if file.endswith('.tsv'):
        temp = pd.read_table('/'.join([root, file]), comment='#', sep = '\t')
        if len(temp.columns) != 5: # using rsem to quantify
            temp = temp.iloc[:,[0, 1, 2, 5]] # gene_id, length, count
            name = file.split('.tsv')[0]
            temp.rename(columns={temp.columns[-1]:name}, inplace= True)
            
            if df.empty:
                df = temp.copy()
            else:
                df = df.merge(temp, on=['transcript_id', 'gene_id', 'length'])

df = df[~df['transcript_id'].str.contains('.PAR_Y')] # remove PAR_Y gene
df = df[df['transcript_id'].str.contains('ENST')] # only keep ENST 

# generate transcript-level matrix for each tissue
for name in tqdm(set(metadata.tissue)):
    col_id = [col for col in metadata[metadata.tissue == name]['id'].tolist() if col in df.columns]
    col_id = ['transcript_id', 'gene_id', 'length'] + col_id
    temp = df[col_id] # tissue-specific isoform transcription
    temp.loc[:,'avgTPM'] = np.round(temp.iloc[:,3:].mean(axis = 1),2)
    temp.loc[:,'gene_id'] = temp['gene_id'].str.replace(r'\.\d+', '', regex=True) # stable gene id
    temp.loc[:,'transcript_id'] = temp['transcript_id'].str.replace(r'\.\d+', '', regex=True) # stable tx id
    temp = temp.merge(gene_type_v29, on = ['gene_id','transcript_id'], how = 'left')
    temp.to_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/total_RNA-seq/tissue/' + name + '.tsv', sep = '\t', header=True, index=False)