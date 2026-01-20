#!/usr/bin/env python
# coding: utf-8

import pandas as pd

meta = pd.read_csv('/gpfs/gibbs/pi/gerstein/yj329/epiPgene/chromatin/preprocessing/bed.exp.tsv', sep ='\t', header = None)
meta.columns = ['id','assay', 'tissue']
meta.loc[:,'name'] = meta.loc[:,'assay'] + '.' + meta.loc[:,'tissue']

root='/gpfs/gibbs/pi/gerstein/yj329/epi/chromatin/preprocessing/bed.files/'

command_1 = 'zcat '
command_2 = 'sort -k1,1 -k2,2n'
command_3 = 'bedtools merge -i stdin -c 1 -o count'

for name in meta['name']:
    acc_ids = [ root + i + '.bed.gz' for i in meta.loc[meta['name']==name, 'id'].to_list()]
    print(command_1 + ' '.join(acc_ids) + ' | ' + command_2 + ' | ' + command_3 + ' > ' + name + '.bed ')
