#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import subprocess
import os,sys


'''
extract pairing information for fastqs
'''

meta_k = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
meta_c = pd.read_csv('../processed_metadata_control.tsv',sep='\t',index_col=0)
meta = pd.concat((meta_k,meta_c),axis=0)
meta['Paired with'] = list(map(lambda x:x.split('/')[2],meta['Paired with']))
meta.rename(columns=lambda x:x.replace(' ','_'),inplace=True)
pair_dict = {}
for row in meta.itertuples():
    file_ = row.Index
    pair_ = row.Paired_with
    order_ = row.Paired_end
    if order_ == 1:
        pair_dict[file_] = pair_
pd.Series(data=pair_dict).to_csv('../pair_info_fastq.txt',sep='\t',header=False)

