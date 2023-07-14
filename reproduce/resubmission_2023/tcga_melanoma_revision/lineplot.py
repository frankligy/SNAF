#!/opt/anaconda3/envs/sctri_env/bin/python3.7

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'



def plot(data_path,ylim,point,name):
    data = pd.read_csv(data_path,index_col=0)
    plt.style.use('dark_background')
    fig,ax = plt.subplots()
    d = data['Distance [µm]'].tolist()
    n = data['DAPI 512(1)'].tolist()
    s = data['GFP 512(1)']
    c = data['Mono(1)']
    ax.plot(d,n,color='#1D00FF',marker='',linestyle='-',label='DAPI')
    ax.plot(d,s,c='#0DF205',marker='',linestyle='-',label='mNeonGreen')
    ax.plot(d,c,c='#F20505',marker='',linestyle='-',label='Cy5')
    ax.set_ylim(ylim)
    ax.set_xlabel('[{}m]'.format('\u03BC'))
    ax.set_ylabel('Intensity')
    ax.axvline(x=point,linestyle='--',c='white')
    ax.set_title(name)
    plt.savefig(os.path.join(os.path.dirname(data_path),'regenerate_lineplot.pdf'),bbox_inches='tight')
    plt.close()

# # SLC45A2, ref
# plot('/Users/ligk2e/Dropbox/to_frank/SLC45A2/ref/Book1_SLC45A2_W1_z12_scale15_data.csv',ylim=[0,3800],point=19.5,name='SLC45A2_ref')

# # SLC45A2, novel
# plot('/Users/ligk2e/Dropbox/to_frank/SLC45A2/novel/Book1_novel_w1_z20_scale15_data.csv',ylim=[0,4100],point=16.8,name='SLC45A2_novel')

# # SIRPA, ref
# plot('/Users/ligk2e/Dropbox/to_frank/SIRPA/ref/Book1_ref_w4_use_new_scale10_data.csv',ylim=[0,2100],point=18.1,name='SIRPA_ref')

# # SIRPA, novel
# plot('/Users/ligk2e/Dropbox/to_frank/SIRPA/novel/Book1_novel_w1_z31_scale1_data.csv',ylim=[0,9500],point=21.3,name='SIRPA_novel')

# # SEMA6A, ref
# plot('/Users/ligk2e/Dropbox/to_frank/SEMA6A/ref/Book1_ref_w1_z1_scale15_data.csv',ylim=[0,490],point=1.9,name='SEMA6A_ref')

# # SEMA6A, novel
# plot('/Users/ligk2e/Dropbox/to_frank/SEMA6A/novel/Book1_novel_w2_z21_scale10.csv',ylim=[0,2400],point=17.8,name='SEMA6A_novel')


# # MET, ref
# plot('/Users/ligk2e/Dropbox/to_frank/MET/ref/Book1_ref_w1_z20_scale15_data.csv',ylim=[0,3900],point=25.8,name='MET_ref')


# MET, novel
plot('/Users/ligk2e/Dropbox/to_frank/MET/novel/Book1_novel_w3_z15_scale15_data.csv',ylim=[0,1350],point=15,name='MET_novel')

