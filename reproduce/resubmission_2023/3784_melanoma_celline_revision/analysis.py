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


# # get neoantigens whose junction are present in this cell line
# tcga = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/T_candidates/T_antigen_candidates_all.txt',sep='\t',index_col=0)
# hla = ['HLA-A*03:01','HLA-A*01:01','HLA-B*07:02','HLA-C*07:02']
# tcga = tcga.loc[tcga['hla'].isin(hla),:]
# df = snaf.remove_trailing_coord('altanalyze_output/ExpressionInput/counts.original.txt',sep='\t')
# uid_cl = set(df.index)
# tcga_subset = tcga.loc[tcga['uid'].isin(uid_cl),:]

# # write and generate fasta
# tcga_subset.to_csv('tcga_melanoma_present.txt',sep='\t')
# with open('search_db.fasta','w') as f:
#     for aa,sub_df in tcga_subset.groupby(by='peptide'):
#         f.write('>' + ','.join(list(set(sub_df['uid'].tolist()))) + '\n')
#         f.write('{}\n'.format(aa))

# configure maxquant
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/melanoma_ms/3784/search_db.fasta']
inputs = [os.path.join('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/melanoma_ms/3784/raw','NCI_3784Mel_Epitope_Rep{}.raw'.format(i+1)) for i in range(5)]
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/melanoma_ms/3784/raw'
snaf.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                   outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)











