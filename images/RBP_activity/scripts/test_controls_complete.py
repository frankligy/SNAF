#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import sys,os

# read raw metadata, test whether all the control samples are in control metadata
kd_metadata_o = pd.read_csv('../metadata_KD.tsv',sep='\t',index_col=0)  # original
kd_metadata_n = pd.read_csv('../metadata_supplement_HNRNPK_FUS_UPF1.tsv',sep='\t',index_col=0)  # new (3 SF that were misclassfied as TF)
kd_metadata = pd.concat((kd_metadata_o,kd_metadata_n),axis=0)
control_metadata = pd.read_csv('../metadata_control.tsv',sep='\t',index_col=0)

control_id = []
for item in kd_metadata['Controlled by']:
    for control in item.split(','):
        control_id.append(control.split('/')[2])
control_id = list(set(control_id))
all_control_id = control_metadata.index.values.tolist()

yes,no = [], []
for item in control_id:
    if item in all_control_id:
        yes.append(item)
    else:
        no.append(item)


# we found no list is not empty, and the root cause is four malformated samples, so delete them from kd metadata
invalid = ['ENCFF588DJP','ENCFF990GRS','ENCFF662AFA','ENCFF697JUO']
processed_kd_metadata = kd_metadata.loc[~(kd_metadata.index.isin(invalid)),:]
processed_kd_metadata.to_csv('../processed_metadata_KD.tsv',sep='\t')

# checked again, now no list is empty
kd_metadata = processed_kd_metadata
control_metadata = pd.read_csv('../metadata_control.tsv',sep='\t',index_col=0)

control_id = []
for item in kd_metadata['Controlled by']:
    for control in item.split(','):
        control_id.append(control.split('/')[2])
control_id = list(set(control_id))
all_control_id = control_metadata.index.values.tolist()

yes,no = [], []
for item in control_id:
    if item in all_control_id:
        yes.append(item)
    else:
        no.append(item)


# only keep the control samples that will be used.
valid = yes
processed_control_metadata = control_metadata.loc[control_metadata.index.isin(valid),:]
processed_control_metadata.to_csv('../processed_metadata_control.tsv',sep='\t')



