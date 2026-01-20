#!/usr/bin/env python
# coding: utf-8
import pandas as pd

meta = pd.read_csv('/gpfs/gibbs/pi/gerstein/yj329/epi/DNAme/preprocessing/bed.exp.tsv', sep ='\t', header = None)
meta.columns = ['id','assay', 'tissue']
meta.loc[:,'name'] = meta.loc[:,'assay'] + '.' + meta.loc[:,'tissue']

root='/gpfs/gibbs/pi/gerstein/yj329/epi/DNAme/preprocessing/bed.files/'
command_1 = 'zcat '
command_2 = 'sort -k1,1 -k2,2n'
command_3 = 'bedtools merge -c 5 -o mean -i - '
for name in meta['name']:
    acc_ids = [ root + i + '.bed.gz' for i in meta.loc[meta['name']==name, 'id'].to_list()]
    print(command_1 + ' '.join(acc_ids) + ' | ' + command_2 + ' | ' + "awk 'NF == 5 && $NF !~ /NaN/ && $NF >= 0 && $NF <= 1' - | " +
          command_3 + '> tissue/' + name + '.bed')
