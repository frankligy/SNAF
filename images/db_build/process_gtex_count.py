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


# argument
data_folder = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision'

# GTEx metadata
sra_table = pd.read_csv(os.path.join(data_folder,'GTEx_2702','GTEx_SRARunTable.txt'),sep='\t')
sra_table = sra_table.loc[:,['Run','body_site']]   # 24455 non-redundant SRR

# # GTEx count
# sra_data = pd.read_csv(os.path.join(data_folder,'GTEx_2702','counts.original.txt'),sep='\t',index_col=0)  
# sra_data.columns.name = None
# sra_data.index.name = None
# sra_data = sra_data.loc[np.logical_not(sra_data.index.duplicated()).tolist(),:]
# sra_data.columns = [item.split('_')[0] for item in sra_data.columns]  
# adata = ad.AnnData(X=csr_matrix(sra_data.values),var=pd.DataFrame(index=sra_data.columns),obs=pd.DataFrame(index=sra_data.index))  
# adata.write(os.path.join(data_folder,'GTEx_2702','gtex_count_raw.h5ad'))
adata = ad.read_h5ad(os.path.join(data_folder,'GTEx_2702','gtex_count_raw.h5ad'))

# only need the overlapping samples 
common = list(set(sra_table['Run'].values).intersection(set(adata.var_names)))  # 2702
sra_table = sra_table.loc[sra_table['Run'].isin(common),:]
srr_to_tissue = sra_table.set_index('Run').squeeze().to_dict()
adata.var['tissue'] = adata.var_names.map(srr_to_tissue).values   
remove_list = ['Cells - Leukemia cell line (CML)','Cells - EBV-transformed lymphocytes','Cells - Cultured fibroblasts']
adata = adata[:,~adata.var['tissue'].isin(remove_list)]  
adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
adata.obs['std'] = adata.X.toarray().std(axis=1)
total_count = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
adata.var['total_count'] = total_count
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
adata_final.write(os.path.join(data_folder,'GTEx_2702','GTEx_junction_counts.h5ad'))
sys.exit('stop')

# inspect
adata = ad.read_h5ad(os.path.join(data_folder,'GTEx_2702','GTEx_junction_counts.h5ad'))
adata.var['tissue'].value_counts().plot(kind='barh',fontsize=2)
plt.savefig(os.path.join(data_folder,'GTEx_2702','tissue_types.pdf'),bbox_inches='tight')
plt.close()













