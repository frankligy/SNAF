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

# GTEx metadata
sra_table = pd.read_csv(os.path.join(data_folder,'scripts','GTEx_SRARunTable.txt'),sep='\t')
sra_table = sra_table.loc[:,['Run','body_site']]   # 24455 non-redundant SRR

# GTEx count
# sra_data = pd.read_csv(os.path.join(data_folder,'GTEx_2702','counts.original.txt'),sep='\t',index_col=0)  
# sra_data.columns.name = None
# sra_data.index.name = None
# sra_data = sra_data.loc[np.logical_not(sra_data.index.duplicated()).tolist(),:]
# sra_data.columns = [item.split('_')[0] for item in sra_data.columns]  
# adata = ad.AnnData(X=csr_matrix(sra_data.values),var=pd.DataFrame(index=sra_data.columns),obs=pd.DataFrame(index=sra_data.index))  
# adata.write(os.path.join(data_folder,'GTEx_2702','gtex_count_raw.h5ad'))
adata = ad.read_h5ad(os.path.join(data_folder,'GTEx_2702','gtex_count_raw.h5ad'))

# only need the overlapping samples 
common = list(set(sra_table['Run'].values).intersection(set(adata.var_names)))
print(len(common))
sra_table = sra_table.loc[sra_table['Run'].isin(common),:]
srr_to_tissue = sra_table.set_index('Run').squeeze().to_dict()
adata.var['tissue'] = adata.var_names.map(srr_to_tissue).values   
remove_list = ['Cells - Leukemia cell line (CML)','Cells - EBV-transformed lymphocytes','Cells - Cultured fibroblasts']
adata = adata[:,~adata.var['tissue'].isin(remove_list)]  
adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
adata.obs['std'] = adata.X.toarray().std(axis=1)
total_count = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
adata.var['total_count'] = total_count
adata.write(os.path.join(data_folder,'GTEx_2702','GTEx_junction_counts.h5ad'))













