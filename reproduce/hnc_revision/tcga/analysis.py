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

# get reduced count matrix and setting up global variable
samples_df = pd.read_csv('sample_hla.txt',sep='\t',index_col=0)
samples = [item.split('_1')[0] + '.bed' for item in samples_df.index.tolist()]
samples_df.index = samples
df = snaf.get_reduced_junction_matrix(pc='counts.TCGA-HNSCC.txt',pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',samples=samples)
sample_to_hla = samples_df.squeeze().to_dict()
hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]

# run SNAF
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
surface.initialize(db_dir=db_dir)

# # running T antigen program
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=30,add_control=add_control,outdir='result',filter_mode='maxmin')  
# jcmq.run(hlas=hlas,strict=False,outdir='./result')
# snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')


# running B antigen program
# membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df,add_control=add_control,not_in_db=False,outdir='result/surface',filter_mode='maxmin')
# surface.run(uids=membrane_tuples,outdir='result/surface',prediction_mode='long_read',n_stride=2,
#             gtf='/data/salomonis2/LabFiles/Frank-Li/refactor/data/2021UHRRIsoSeq_SQANTI3_filtered.gtf',
#             tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',serialize=True)
surface.generate_full_results(outdir='result/surface',freq_path='result/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt',mode='long_read',validation_gtf=None)

# short-read mode
membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df,add_control=add_control,not_in_db=False,outdir='result/surface',filter_mode='maxmin')
surface.run(uids=membrane_tuples,outdir='result/surface',prediction_mode='short_read',n_stride=2,
            gtf=None,
            tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',serialize=True)
surface.generate_full_results(outdir='result/surface',freq_path='result/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt',mode='short_read',validation_gtf='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/SQANTI-all/collapse_isoforms_classification.filtered_lite.gtf')
