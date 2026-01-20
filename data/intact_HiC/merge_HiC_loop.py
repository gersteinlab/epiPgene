#!/usr/bin/env python

import pandas as pd
meta = pd.read_csv('/gpfs/gibbs/pi/gerstein/yj329/epi/HiC/preprocessing/bedpe.exp.tsv', sep ='\t', header = None)
meta.columns = ['id','assay', 'tissue']

root='/gpfs/gibbs/pi/gerstein/yj329/epi/HiC/preprocessing/loop.files/'
for name in set(meta.tissue):
    acc_ids = [ root + i + '.bedpe.gz' for i in meta.loc[meta['tissue']==name, 'id'].to_list()]
    print('zcat ' + ' '.join(acc_ids) + " | grep -v '#' | cut -f1-24 | sort -k1,1 -k2,2n > " + name + '.bedpe')