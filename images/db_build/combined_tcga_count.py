#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt
import anndata as ad
from scipy.sparse import csr_matrix
from tqdm import tqdm
import pickle


cancers = ['BLCA','BRCA','COAD','ESCA','GBM','HNSCC','KICH','KIRC','KIRP','LIHC','Lung','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC']
root = '/data/salomonis2/NCI-R01/control-counts_TCGA'


# events = set()
# for c in tqdm(cancers):
#     path = os.path.join(root,'counts.TCGA-{}-controls.txt'.format(c))
#     df = pd.read_csv(path,sep='\t',index_col=0) 
#     events = events.union(set(df.index))
# events = list(events)   # 6552446
# with open('tcga_events_all.p','wb') as f:
#     pickle.dump(events,f)


df_list = []
id_list = []
for c in tqdm(cancers):
    path = os.path.join(root,'counts.TCGA-{}-controls.txt'.format(c))
    df = pd.read_csv(path,sep='\t',index_col=0)
    if df.shape[1] > 0:
        df_list.append(df)
        id_list.extend([c]*df.shape[1])
combined_df = pd.concat(df_list,axis=1,join='outer').fillna(0)
adata = ad.AnnData(X=csr_matrix(combined_df.values),obs=pd.DataFrame(index=combined_df.index),var=pd.DataFrame(index=combined_df.columns))
adata.var['tissue'] = id_list
adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
adata.obs['std'] = adata.X.toarray().std(axis=1)
total_count = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
adata.var['total_count'] = total_count
print(adata)
adata.write('../TCGA_matched_control/tcga_matched_control_junction_count.h5ad')
adata.var['tissue'].value_counts().plot(kind='barh')
plt.savefig('../TCGA_matched_control/tissue_types.pdf',bbox_inches='tight')
plt.close()














