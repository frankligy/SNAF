#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os,sys
sys.path.insert(0,'/data/salomonis2/software')
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import subprocess
import anndata as ad

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# get reduced count matrix and setting up global variable
df = snaf.get_reduced_junction_matrix(pc='counts.original.txt',pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')  
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
surface.initialize(db_dir=db_dir)

# # running T antigen program
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=16,add_control=add_control,outdir='result',filter_mode='maxmin')
# sample_to_hla = pd.read_csv('sample_hla_new.txt',sep='\t',index_col=0)['hla'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq.run(hlas=hlas,strict=False,outdir='./result')
# snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')

# look at a few t antigen
bam_path_list = ['/data/salomonis-archive/FASTQs/NCI-R01/Wise-Draper/Frank_New_analysis/bams/RNA_DE-0404832_S26_secondAligned.sortedByCoord.out.bam',
                 '/data/salomonis-archive/FASTQs/NCI-R01/Wise-Draper/Frank_New_analysis/bams/RNA_TN20-158485_S91_secondAligned.sortedByCoord.out.bam',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-HNSCC/TCGA-CV-A45Y-01A.bam',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-HNSCC/TCGA-P3-A6T8-01A.bam']
bai_path_list = ['/data/salomonis-archive/FASTQs/NCI-R01/Wise-Draper/Frank_New_analysis/bams/RNA_DE-0404832_S26_secondAligned.sortedByCoord.out.bam.bai',
                 '/data/salomonis-archive/FASTQs/NCI-R01/Wise-Draper/Frank_New_analysis/bams/RNA_TN20-158485_S91_secondAligned.sortedByCoord.out.bam.bai',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-HNSCC/TCGA-CV-A45Y-01A.bai',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-HNSCC/TCGA-P3-A6T8-01A.bai']
outdir = '../Frank_inspection/ggsashimi_run_folder'
sif_anno_path = '/data/salomonis2/software/ggsashimi'
bam_contig_rename = [True,True,False,False]

'''
tested uid:
1. ENSG00000148219:E63.4_116840613-E73.1
2. ENSG00000117122:U0.1_16981514-I2.2
3. ENSG00000166415:E23.3-E24.1
4. ENSG00000170298:U0.1_20468392-E2.1
5. ENSG00000187689:E2.1-E3.2         switch tumor tcga sample to TCGA-CV-A45Y-01A and TCGA-P3-A6T8-01A, before we use TCGA-CN-4739-01A and TCGA-CN-6988-01A
6. ENSG00000136231:E3.1-E6.1
7. ENSG00000162892:E2.2-E3.2
8. ENSG00000136231:E6.1-E14.1
'''

uid = 'ENSG00000136231:E6.1-E14.1'
snaf.gtex_visual_combine_plotly(uid=uid,outdir='../Frank_inspection',norm=False,tumor=df)
snaf.gtex_visual_combine_plotly(uid=uid,outdir='../Frank_inspection',norm=True,tumor=df)
snaf.prepare_sashimi_plot(bam_path_list,bai_path_list,outdir,sif_anno_path,bam_contig_rename,query_region='chr7:23361541-23418976',skip_copy=True, min_junction=20)
sys.exit('stop')

'''B cell neoantigen'''
# long-read mode
# membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df,add_control=add_control,not_in_db=False,outdir='result/surface',filter_mode='maxmin')
# surface.run(uids=membrane_tuples,outdir='result/surface',prediction_mode='long_read',n_stride=2,
#             gtf='/data/salomonis2/LabFiles/Frank-Li/refactor/data/2021UHRRIsoSeq_SQANTI3_filtered.gtf',
#             tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',serialize=True)
# surface.generate_full_results(outdir='result/surface',freq_path='result/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt',mode='long_read',validation_gtf=None)

# short-read mode
# membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df,add_control=add_control,not_in_db=False,outdir='result/surface',filter_mode='maxmin')
# surface.run(uids=membrane_tuples,outdir='result/surface',prediction_mode='short_read',n_stride=2,
#             gtf=None,
#             tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',serialize=True)
# surface.generate_full_results(outdir='result/surface',freq_path='result/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt',mode='short_read',validation_gtf='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/SQANTI-all/collapse_isoforms_classification.filtered_lite.gtf')




