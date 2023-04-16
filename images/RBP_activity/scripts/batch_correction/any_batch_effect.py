#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os
import sys
import subprocess
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from sklearn.preprocessing import MinMaxScaler
import umap
import matplotlib.patches as mpatches
sys.setrecursionlimit(1000000)

from pca_umap import *

# some input needed
query = 'shRNA_HepG2'
event_annotation_file = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/limma_bcAltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt'
batch_file = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/scripts/batch_correction/batch_shRNA_HepG2.txt'
group_file = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/scripts/batch_correction/groups_exp_control.txt'
plot_outdir = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/batch_effect/shRNA_HepG2'
intra_correlation = 0.8

if not os.path.exists(plot_outdir):
    os.mkdir(plot_outdir)
batch = pd.read_csv(batch_file,sep='\t',index_col=0)
batch_dict = batch.squeeze().to_dict()
group = pd.read_csv(group_file,sep='\t',header=None)
group_dict = group.set_index(keys=0)[2].to_dict()
df = pd.read_csv(event_annotation_file,sep='\t')
df_uid = df['UID'].to_frame()
df_samples = df.loc[:,batch.index]
df = pd.concat([df_uid,df_samples],axis=1).set_index('UID')  
fill_as_zero = df.fillna(value=0)
def median_impute(x):
    med = np.ma.median(np.ma.masked_invalid(x.values))
    result = np.nan_to_num(x.values,nan=med)
    return result
fill_as_median = df.apply(func=median_impute,axis=1,result_type='expand')


adata = ad.AnnData(X=fill_as_median.values,var=pd.DataFrame(index=df.columns.values),obs=pd.DataFrame(index=df.index.values))
adata.layers['compute_prop'] = fill_as_zero.values
adata.var['batch'] = adata.var_names.map(batch_dict).values
adata.var['is_control'] = adata.var_names.map(group_dict).values
adata.obs['gene'] = [item[0] for item in adata.obs_names.str.split(':')]
proportion = np.count_nonzero(adata.layers['compute_prop'],axis=1) / adata.layers['compute_prop'].shape[1]
adata.obs['proportion'] = proportion
adata.obs['cv'] = np.std(adata.X,axis=1)/np.mean(adata.X,axis=1)

## intra-correlation > cutoff should be removed to avoid collinearity
adata.obs['order'] = np.arange(adata.obs.shape[0])
new_obs_chunks = []
for gene, grouped_df in adata.obs.groupby(by='gene'):
    mat = adata[grouped_df.index,:].X
    try:
        u_tri_mat = np.triu(np.absolute(np.corrcoef(mat)),k=1)
    except TypeError: 
        # TypeError: tri() missing 1 required positional argument: 'N'
        # that is because only 1 event in the gene, np.triu(1) will raise error
        grouped_df['to_drop_{}'.format(intra_correlation)] = [False]
        new_obs_chunks.append(grouped_df)
        continue
    to_drop = []
    for j in range(u_tri_mat.shape[1]):
        if np.any(u_tri_mat[:,j] > intra_correlation):
            to_drop.append(True)
        else:
            to_drop.append(False)
    grouped_df['to_drop_{}'.format(intra_correlation)] = to_drop
    new_obs_chunks.append(grouped_df)
adata.obs = pd.concat(new_obs_chunks).sort_values(by='order')

draw_pca_umap(adata,'shRNA_K562',outdir=plot_outdir)
for column in ['batch','is_control']:
    draw_all_pca_umap(adata,'shRNA_K562',column=column,outdir=plot_outdir,layer=None)





