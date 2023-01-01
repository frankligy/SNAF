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
df = snaf.get_reduced_junction_matrix(pc='altanalyze_output/ExpressionInput/counts.original.txt',
                                      pea='altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')

# run SNAF
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
surface.initialize(db_dir=db_dir)

'''T cell neoantigen'''
jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=16,add_control=add_control,outdir='result',filter_mode='maxmin')
hlas = [['HLA-A*02:01', 'HLA-B*15:11', 'HLA-C*03:03'],['HLA-A*02:01', 'HLA-B*15:11', 'HLA-C*03:03'],['HLA-A*02:01', 'HLA-B*15:11', 'HLA-C*03:03']]
jcmq.run(hlas=hlas,outdir='./result')
snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')
sys.exit('stop')

'''B cell neoantigen'''
# long-read
membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df,add_control=add_control,not_in_db=False,outdir='result_new/surface',filter_mode='maxmin')
surface.run(uids=membrane_tuples,outdir='result_new/surface',prediction_mode='long_read',n_stride=2,
            gtf='/data/salomonis2/LabFiles/Frank-Li/refactor/data/2021UHRRIsoSeq_SQANTI3_filtered.gtf',
            tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',serialize=True)
surface.generate_full_results(outdir='result_new/surface',freq_path='result_new/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt',mode='long_read',validation_gtf=None)

# short-read
membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df,add_control=add_control,not_in_db=False,outdir='result_new/surface',filter_mode='maxmin')
surface.run(uids=membrane_tuples,outdir='result_new/surface',prediction_mode='short_read',n_stride=2,
            gtf=None,
            tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',serialize=True)
surface.generate_full_results(outdir='result_new/surface',freq_path='result_new/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt',mode='short_read',validation_gtf='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/SQANTI-all/collapse_isoforms_classification.filtered_lite.gtf')

# surface.run_dash_B_antigen(pkl='result/surface/surface_antigen_lr.p',candidates='result/surface/candidates_3_lr_None_True.txt',prediction_mode='long_read',
#                            python_executable='/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7')



sys.exit('stop')



### MS validation
# df = pd.read_csv('result/frequency_stage3_verbosity1_uid_gene_symbol.txt',sep='\t',index_col=0)
# df_common = df.loc[df['n_sample']>282,:]
# df_unique = df.loc[df['n_sample']==1,:]
# with open('result/MS_validation_common.fasta','w') as f:
#     for row in df_common.itertuples():
#         peptide, uid = row.Index.split(',')
#         f.write('>{}\n{}\n'.format(uid,peptide))
# with open('result/MS_validation_unique.fasta','w') as f:
#     for row in df_unique.itertuples():
#         peptide, uid = row.Index.split(',')
#         f.write('>{}\n{}\n'.format(uid,peptide))
# for f in ['MS_validation_common','MS_validation_unique']:
#     snaf.proteomics.remove_redundant('./result/{}.fasta'.format(f),'./result/{}_unique.fasta'.format(f))
#     snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot_9_10_mers_unique.fasta',
#                                       fa2_path='./result/{}_unique.fasta'.format(f),outdir='./result',write_unique2=True,prefix='{}_'.format(f))

db_common = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/fasta/MS_validation_common_unique2.fasta']
db_unique = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/fasta/MS_validation_unique_unique2.fasta']

'''
common
'''
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# for p in all_patients:
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common/{}'.format(p))
#     all_samples = subprocess.run('for file in *.raw; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     inputs = []
#     pwd = os.getcwd()
#     for s in all_samples:
#         inputs.append(os.path.join(pwd,s))
#     snaf.proteomics.set_maxquant_configuration(dbs=db_common,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                                outdir=pwd,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)


