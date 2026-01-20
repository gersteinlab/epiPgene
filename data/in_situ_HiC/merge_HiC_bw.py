#!/usr/bin/env python

import pandas as pd
meta = pd.read_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/in_situ_HiC/preprocessing/bigwig.exp.tsv', sep ='\t', header = None)
meta.columns = ['id','assay', 'tissue']

root='/gpfs/gibbs/pi/gerstein/yj329/epiPgene/in_situ_HiC/preprocessing/bigwig.files/'
command_1 = 'wiggletools mean '
command_2 = 'wigToBigWig -clip '
for name in set(meta.tissue):
    acc_ids = [ root + i + '.bigWig' for i in meta.loc[meta['tissue']==name, 'id'].to_list()]
    print(command_1 + ' '.join(acc_ids) + ' > ' + name + '.wig ; ' + 
          command_2 + name + '.wig ' + 'hg38.chrom.sizes ' + name + '.bw ; ' +
         'rm ' + name + '.wig')