#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import subprocess
import os,sys


'''
This script is to compare the overlap between the 299 SF in KD and the SF in binding side
'''

os.chdir('/data/salomonis2/NCI-R01/Signatures/encode_eclip_bed/raw/K562/idr_bed')
bind_k562 = subprocess.run("for file in *.bed; do echo $file; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
bind_k562 = [item.split('.')[0] for item in bind_k562]   # 120
os.chdir('/data/salomonis2/NCI-R01/Signatures/encode_eclip_bed/raw/HepG2/idr_bed')
bind_hepg2 = subprocess.run("for file in *.bed; do echo $file; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
bind_hepg2 = [item.split('.')[0] for item in bind_hepg2]  # 103
common = list(set(bind_k562).intersection(set(bind_hepg2)))  # 73
union = list(set(bind_k562).union(set(bind_hepg2)))   # 150

os.chdir('/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/scripts')
ori_kd_meta = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
supp_kd_meta = pd.read_csv('../metadata_supplement_HNRNPK_FUS_UPF1.tsv',sep='\t',index_col=0)
concat_kd_meta = pd.concat([ori_kd_meta,supp_kd_meta],axis=0)
kd_all = list(set([item.split('-')[0] for item in concat_kd_meta['Experiment target']]))  # 299

overlap = list(set(union).intersection(set(kd_all)))  # 125







