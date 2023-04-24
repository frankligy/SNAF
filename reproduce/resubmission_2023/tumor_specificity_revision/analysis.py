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


# get reduced junction
df = snaf.get_reduced_junction_matrix(pc='counts.TCGA-SKCM.txt',pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')

# run SNAF
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
surface.initialize(db_dir=db_dir)

# # output the combined adata for melanoma
# snaf.gtex.adata.write('combined_normal_count.h5ad')

# plot
uid = 'ENSG00000128422:E5.3-E6.1'
snaf.gtex_visual_subplots(uid=uid,norm=True,outdir='./')
snaf.gtex_visual_per_tissue_count(uid=uid,outdir='./')
sys.exit('stop')

# testing bayesian statistics
'''
ENSG00000105976:E3.1-E4.1
ENSG00000198053:E7.2-E13.1_1915159
ENSG00000164175:E3.2_33963931-E4.2
ENSG00000092421:E22.1-E24.1_116468915
ENSG00000057019:E5.1-I5.1_98867438
ENSG00000152558:I7.1_102419074-E8.1
'''



