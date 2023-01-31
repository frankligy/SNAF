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


# # load GTEx junction-level gene count
# path = '/data/salomonis2/LabFiles/Frank-Li/refactor/GTEx_reprocess/altanalyze_output/ExpressionInput/counts.original-steady-state.txt'
# sra_data = pd.read_csv(path,sep='\t',index_col=0)
# sra_data.columns.name = None
# sra_data.index.name = None
# sra_data = sra_data.loc[np.logical_not(sra_data.index.duplicated()).tolist(),:]
# sra_data.columns = [item.split('_')[0] for item in sra_data.columns]  
# adata = ad.AnnData(X=csr_matrix(sra_data.values),var=pd.DataFrame(index=sra_data.columns),obs=pd.DataFrame(index=sra_data.index))  
# data_folder = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision'
# sra_table = pd.read_csv(os.path.join(data_folder,'GTEx_2702','GTEx_SRARunTable.txt'),sep='\t')
# sra_table = sra_table.loc[:,['Run','body_site']]   # 24455 non-redundant SRR
# common = list(set(sra_table['Run'].values).intersection(set(adata.var_names)))  # 2702
# sra_table = sra_table.loc[sra_table['Run'].isin(common),:]
# srr_to_tissue = sra_table.set_index('Run').squeeze().to_dict()
# adata.var['tissue'] = adata.var_names.map(srr_to_tissue).values   
# remove_list = ['Cells - Leukemia cell line (CML)','Cells - EBV-transformed lymphocytes','Cells - Cultured fibroblasts']
# adata = adata[:,~adata.var['tissue'].isin(remove_list)]  
# adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
# adata.obs['std'] = adata.X.toarray().std(axis=1)
# total_count = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
# adata.var['total_count'] = total_count
# adata.write('gtex_count.h5ad') # 33856 × 2629

# load TCGA matched controls
# cancers = ['BLCA','BRCA','COAD','ESCA','GBM','HNSCC','KICH','KIRC','KIRP','LIHC','Lung','OV','PAAD','PCPG','PRAD','READ','SARC','SKCM','STAD','TGCT','THCA','THYM','UCEC']
# root = '/data/salomonis2/NCI-R01/ONCObrowser/GeneExpression/Counts'
# df_list = []
# id_list = []
# for c in tqdm(cancers):
#     path = os.path.join(root,'counts.TCGA-{}-steady-state_controls.txt'.format(c))
#     df = pd.read_csv(path,sep='\t',index_col=0)
#     if df.shape[1] > 0:
#         df_list.append(df)
#         id_list.extend([c]*df.shape[1])
# combined_df = pd.concat(df_list,axis=1,join='outer').fillna(0)
# adata = ad.AnnData(X=csr_matrix(combined_df.values),obs=pd.DataFrame(index=combined_df.index),var=pd.DataFrame(index=combined_df.columns))
# adata.var['tissue'] = id_list
# adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
# adata.obs['std'] = adata.X.toarray().std(axis=1)
# total_count = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
# adata.var['total_count'] = total_count    # 36911 × 702
# adata.write('tcga_matched_control_count.h5ad')

# # load gtex skin control
# path = '/data/salomonis-archive/BAMs/NCI-R01/GTEx/newSkin_100922/Frank/altanalyze_output/ExpressionInput/counts.original-steady-state.txt'
# sra_data = pd.read_csv(path,sep='\t',index_col=0)
# sra_data.columns.name = None
# sra_data.index.name = None
# sra_data = sra_data.loc[np.logical_not(sra_data.index.duplicated()).tolist(),:]
# sra_data.columns = [item.split('_')[0] for item in sra_data.columns]  
# adata = ad.AnnData(X=csr_matrix(sra_data.values),var=pd.DataFrame(index=sra_data.columns),obs=pd.DataFrame(index=sra_data.index))  
# adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
# adata.obs['std'] = adata.X.toarray().std(axis=1)
# total_count = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
# adata.var['total_count'] = total_count  # 24962 × 313
# adata.write('gtex_new_skin_count.h5ad')

def gtex_configuration(df,gtex_db,add_control=None):
    global adata_gtex
    global adata
    global t_min
    global n_max
    global normal_cutoff
    global tumor_cutoff
    global normal_prevalance_cutoff
    global tumor_prevalance_cutoff
    tested_junctions = set(df.index)
    adata = ad.read_h5ad(gtex_db)
    adata = adata[np.logical_not(adata.obs_names.duplicated()),:] 
    adata = adata[list(set(adata.obs_names).intersection(tested_junctions)),:]  
    print('Current loaded gtex cohort with shape {}'.format(adata.shape))
    tissue_dict = adata.var['tissue'].to_dict()
    adata_gtex = adata   # already has mean and tissue variables
    if add_control is not None:
        for id_, control in add_control.items():
            if isinstance(control,pd.DataFrame):
                assert len(set(control.columns).intersection(tissue_dict.keys())) == 0  # sample id can not be ambiguous
                control = control.loc[np.logical_not(control.index.duplicated()),:]
                control = control.loc[list(set(control.index).intersection(tested_junctions)),:]
                print('Adding cohort {} with shape {} to the database'.format(id_,control.shape))
                tissue_dict_right = {k:id_ for k in control.columns}
                tissue_dict.update(tissue_dict_right)
                df_left = adata.to_df()
                df_right = control
                df_combine = df_left.join(other=df_right,how='outer').fillna(0)
                adata = ad.AnnData(X=csr_matrix(df_combine.values),obs=pd.DataFrame(index=df_combine.index),var=pd.DataFrame(index=df_combine.columns))

            elif isinstance(control,ad.AnnData):
                assert len(set(control.var_names).intersection(tissue_dict.keys())) == 0
                control = control[np.logical_not(control.obs_names.duplicated()),:]
                control = control[list(set(control.obs_names).intersection(tested_junctions)),:]
                print('Adding cohort {} with shape {} to the database'.format(id_,control.shape))
                if 'tissue' in control.var.columns:   # if tissue is in var columns, it will be used 
                    tissue_dict_right = control.var['tissue'].to_dict()
                else:
                    tissue_dict_right = {k:id_ for k in control.var_names}
                tissue_dict.update(tissue_dict_right)
                df_left = adata.to_df()
                df_right = control.to_df()
                df_combine = df_left.join(other=df_right,how='outer').fillna(0)
                adata = ad.AnnData(X=csr_matrix(df_combine.values),obs=pd.DataFrame(index=df_combine.index),var=pd.DataFrame(index=df_combine.columns))
            
            else:
                raise Exception('control must be either in dataframe or anndata format')

            print('now the shape of control db is {}'.format(adata.shape))
            adata.var['tissue'] = adata.var_names.map(tissue_dict).values
            adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
            total_count = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
            adata.var['total_count'] = total_count
            


    return adata

tcga_ctrl_db = ad.read_h5ad('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ts/gene/tcga_matched_control_count.h5ad')
gtex_skin_ctrl_db = ad.read_h5ad('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ts/gene/gtex_new_skin_count.h5ad')
add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}
gtex_db = ad.read_h5ad('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ts/gene/gtex_count.h5ad')
union = list(set(gtex_db.obs_names).union(set(tcga_ctrl_db.obs_names)).union(set(gtex_skin_ctrl_db.obs_names))) # 38053
fake_df = pd.DataFrame(index=union)
adata = gtex_configuration(df=fake_df,
                           gtex_db='/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ts/gene/gtex_count.h5ad',
                           add_control=add_control)
adata.write('combined_gene_count.h5ad')  # 38053 × 3644
                         


