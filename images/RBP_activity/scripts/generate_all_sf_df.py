#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import subprocess
import os,sys


'''
This script is to generate a txt file for every SF, I want to know how their control samples were connected with each other.
'''

ori_kd_meta = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
supp_kd_meta = pd.read_csv('../metadata_supplement_HNRNPK_FUS_UPF1.tsv',sep='\t',index_col=0)
concat_kd_meta = pd.concat([ori_kd_meta,supp_kd_meta],axis=0)
need_kd_meta = concat_kd_meta.loc[:,['Assay','Biosample term name','Experiment target','Biological replicate(s)','Technical replicate(s)','Paired end','Controlled by']]
for sf,sub_df in need_kd_meta.groupby(by='Experiment target'):
    sub_df.to_csv('../each_sf_df/{}_sub_df.txt'.format(sf),sep='\t')




