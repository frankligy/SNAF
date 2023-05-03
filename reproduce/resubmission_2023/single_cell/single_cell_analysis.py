#!/data/salomonis2/LabFiles/Frank-Li/scTriangulate/spatial_slide/sctri_spatial_env/bin/python3.7

import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import subprocess
import os,sys
sys.path.insert(0,'/data/salomonis2/software')
from sctriangulate import *
from sctriangulate.preprocessing import *
from sctriangulate.colors import bg_greyed_cmap


mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# adata = small_txt_to_adata(int_file='normed_expression_matrix.txt',gene_is_index=True,sep='\t')
# adata = scanpy_recipe(adata,is_log=True,resolutions=[1],modality='rna',pca_n_comps=50,n_top_genes=3000)

adata = sc.read('adata_after_scanpy_recipe_rna_1_umap_True.h5ad')
add_annotations(adata,inputs='celltype.txt',cols_input=['ct','tumor'],index_col=0,cols_output=['ct','tumor'],kind='disk')
# umap_dual_view_save(adata,cols=['ct','tumor'])
pd.DataFrame(data=adata.obsm['X_umap'],index=adata.obs_names,columns=['umap_x','umap_y']).to_csv('umap_coord.txt',sep='\t')
sys.exit('stop')

# map psi value of certain splicing events
lookup = pd.read_csv('lookup.txt',sep='\t',index_col=0).set_index(keys='barcode')
lookup = lookup['fastq'].to_dict()
reverse_lookup = {v:k for k,v in lookup.items()} 

# tumor spefiic
# splicing_events = [
#     'MYL6:ENSG00000092841:E3.12-E3.15|ENSG00000092841:E3.12-E4.1',
#     'PSAP:ENSG00000197746:E8.1-E9.3|ENSG00000197746:E7.1-E9.3'
# ]

# tcell spefiic
# splicing_events = [
#     'PTPRC:ENSG00000081237:E20.1-E30.1|ENSG00000081237:E25.2-E30.1',
#     'HSPD1:ENSG00000144381:E3.3-I3.1|ENSG00000144381:E2.12-E4.1'
# ]
# ee = pd.read_csv('altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t',index_col='UID').iloc[:,10:]
# added_df = ee.loc[splicing_events,:]
# renamed_column = [reverse_lookup[item.split('_secondAligned.sortedByCoord.out.bed')[0]] for item in added_df.columns]
# added_df.columns = renamed_column
# added_df = added_df.T.fillna(0)

# adata = adata[renamed_column,:]

# add_annotations(adata,inputs=added_df,cols_input=splicing_events,index_col=0,cols_output=[item.replace(':','_').replace('|','_') for item in splicing_events],kind='memory')
# for item in [item.replace(':','_').replace('|','_') for item in splicing_events]:
#     adata.obs[item] = adata.obs[item].astype('float32')
#     sc.pl.umap(adata,color=item,vmin=1e-5,cmap=bg_greyed_cmap('viridis'))
#     plt.savefig('umap_{}.pdf'.format(item),bbox_inches='tight')

# visualize all 542 junction
common = []
with open('common_542_junction.txt','r') as f:
    for line in f:
        common.append(line.rstrip('\n'))
ee = pd.read_csv('altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t',index_col='UID').iloc[:,10:]
junc = [':'.join(item.split('|')[0].split(':')[1:]) for item in ee.index]
ee.index = junc
ee = ee.loc[~ee.index.duplicated(),:]
splicing_events = common
added_df = ee.loc[splicing_events,:]
renamed_column = [reverse_lookup[item.split('_secondAligned.sortedByCoord.out.bed')[0]] for item in added_df.columns]
added_df.columns = renamed_column
added_df = added_df.T.fillna(0)

adata = adata[renamed_column,:]

add_annotations(adata,inputs=added_df,cols_input=splicing_events,index_col=0,cols_output=[item.replace(':','_') for item in splicing_events],kind='memory')
for item in [item.replace(':','_') for item in splicing_events]:
    adata.obs[item] = adata.obs[item].astype('float32')
    sc.pl.umap(adata,color=item,vmin=1e-5,cmap=bg_greyed_cmap('viridis'))
    plt.savefig('common_542_junction_umap/umap_{}.pdf'.format(item),bbox_inches='tight')



