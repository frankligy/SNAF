#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os

meta_k = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
url_k = meta_k['File download URL'].tolist()

meta_c = pd.read_csv('../processed_metadata_control.tsv',sep='\t',index_col=0)
url_c = meta_c['File download URL'].tolist()


with open('../fastq.txt','w') as f:
    for item in url_k + url_c:
        f.write('"{}"\n'.format(item))




