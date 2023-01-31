#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
sys.path.insert(0,'/data/salomonis2/software')
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import anndata as ad
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import subprocess

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
peptide = pd.read_csv('inputs.txt',sep='\t',index_col=0).index.tolist()
hlac = pd.read_csv('hlac.txt',sep='\t',index_col=0).index.tolist()
hlac = snaf.snaf.hla_formatting(hlac,'netMHCpan_output','netMHCpan_input')
dfs = []
for pep in peptide:
    df = snaf.binding.run_netMHCpan(software_path=netMHCpan_path,
                                    peptides=[pep],
                                    hlas=hlac,
                                    length=len(pep),
                                    cmd_num=1,
                                    tmp_dir=None,
                                    tmp_name=None)
    dfs.append(df)
final_df = pd.concat(dfs,axis=0)
final_df.to_csv('prediction_netMHCpan.txt',sep='\t')

final_df = pd.read_csv('prediction_netMHCpan.txt',sep='\t',index_col=0)
final_df_valid = final_df.loc[(final_df['mer'] == 9) | (final_df['mer'] == 10),['peptide','hla']]
final_df_valid['hla'] = snaf.snaf.hla_formatting(final_df_valid['hla'].tolist(),'netMHCpan_output','deepimmuno')
final_df_valid = snaf.deepimmuno.run_deepimmuno(final_df_valid)
final_df_valid['HLA'] = snaf.snaf.hla_formatting(final_df_valid['HLA'].tolist(),'deepimmuno','netMHCpan_output')
final_df['uid'] = [p + ',' + h for p,h in zip(final_df['peptide'],final_df['hla'])]
final_df.set_index('uid',inplace=True)
final_df_valid['uid'] = [p + ',' + h for p,h in zip(final_df_valid['peptide'],final_df_valid['HLA'])]
final_df_valid.set_index('uid',inplace=True)
final_df = final_df.join(final_df_valid,how='left',lsuffix='_left',rsuffix='_right')
final_df.to_csv('add_deepimmuno.txt',sep='\t')


