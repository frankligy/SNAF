#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
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



# construct a metadata table
df_mxd4 = pd.read_csv('GSE202700_SraRunTable.csv',index_col=0)
df_hypomethyl = pd.read_csv('GSE197471_SraRunTable.csv',index_col=0)  # 2*50
df_scalp = pd.read_csv('GSE149189_SraRunTable.csv',index_col=0)   # 1*75
df_klf4 = pd.read_csv('GSE111786_SraRunTable.csv',index_col=0)
df_msc = pd.read_csv('GSE102983_SraRunTable.csv',index_col=0)


samples_dict = {}
for srr,gt,tis in zip(df_mxd4.index,df_mxd4['Genotype'],df_mxd4['Tissue']):
    samples_dict[srr] = '_'.join([tis.replace('\, ','_'),gt.replace(' ','_')])
for srr,cl,oxy in zip(df_hypomethyl.index,df_hypomethyl['Cell_Line'],df_hypomethyl['culture_oxygen']):
    samples_dict[srr] = '_'.join([cl,oxy.replace(' ','_').replace('%','')])
for srr,age,tis in zip(df_scalp.index,df_scalp['Age'],df_scalp['Tissue']):
    samples_dict[srr] = '_'.join([tis.replace(' ','_'),age])
for srr in df_klf4.index:
    samples_dict[srr] = 'skin_epidermis_klf4_study'
for srr,ct in zip(df_msc.index,df_msc['Cell_type']):
    samples_dict[srr] = ct.replace(' ','_')
pd.Series(samples_dict,name='label').to_csv('sample_meta.txt',sep='\t')





















