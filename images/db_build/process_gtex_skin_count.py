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


# argument
data_folder = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision'

# GTEx count
sra_data = pd.read_csv(os.path.join(data_folder,'GTEx_skin','counts.original.txt'),sep='\t',index_col=0)  
sra_data.columns.name = None
sra_data.index.name = None
sra_data = sra_data.loc[np.logical_not(sra_data.index.duplicated()).tolist(),:]
sra_data.columns = [item.split('_')[0] for item in sra_data.columns]  
adata = ad.AnnData(X=csr_matrix(sra_data.values),var=pd.DataFrame(index=sra_data.columns),obs=pd.DataFrame(index=sra_data.index))   # 877509 × 313
adata.obs_names = [item.split('=')[0] for item in adata.obs_names]
## for duplicated junction, only keep the one with highest junction counts
junctions_dup = list(set(adata.obs_names[adata.obs_names.duplicated()]))
junctions_not_dup = list(set(adata.obs_names).difference(set(junctions_dup)))
adata_tmp = adata[np.logical_not(adata.obs_names.duplicated()),:].copy()
adata_not_dup = adata_tmp[junctions_not_dup,:].copy()
adata_dup_stack = []
for junc in tqdm(junctions_dup,total=len(junctions_dup)):
    adata_junc = adata[junc,:].copy()
    row_index_with_max_count = np.argmax(adata_junc.obs['mean'].values)
    adata_dup_stack.append(adata_junc[row_index_with_max_count,:].copy())
adata_dup = ad.concat(adata_dup_stack,axis=0,join='outer',merge='first')
adata_final = ad.concat([adata_not_dup,adata_dup],axis=0,join='outer',merge='first')
adata_final.write(os.path.join(data_folder,'GTEx_skin','gtex_skin_count.h5ad'))
















