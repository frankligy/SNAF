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

'''
raw file: /data/salomonis-archive/FASTQs/NCI-R01/Leukemia/Leucegene-Project
'''

# just see if the neoantigens from the 51 common neojunctions are present
list51 = pd.read_csv('attachment.txt',sep='\t',index_col=0,header=None).index.tolist()
tcga = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/T_candidates/T_antigen_candidates_all.txt',sep='\t',index_col=0)
hla = ['HLA-A*03:01','HLA-A*01:01','HLA-B*07:02','HLA-C*07:02']
tcga_subset = tcga.loc[(tcga['hla'].isin(hla)) & (tcga['uid'].isin(list51)),:]
with open('search_db.fasta','w') as f:
    for aa,sub_df in tcga_subset.groupby(by='peptide'):
        f.write('>' + ','.join(list(set(sub_df['uid'].tolist()))) + '\n')
        f.write('{}\n'.format(aa))


# # run snaf on this dataset

# # get reduced junction
# df = snaf.remove_trailing_coord('altanalyze_output/ExpressionInput/counts.original.txt',sep='\t')

# # run SNAF
# netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
# db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
# tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
# gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
# add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

# snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
# surface.initialize(db_dir=db_dir)

# T antigen
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=16,add_control=add_control,outdir='result',filter_mode='maxmin')
# hla = ['HLA-A*03:01','HLA-A*01:01','HLA-B*07:02','HLA-B*07:02','HLA-C*07:02','HLA-C*07:02']
# hlas = [hla,hla]
# jcmq.run(hlas=hlas,outdir='./result')
# snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')
# jcmq = snaf.JunctionCountMatrixQuery.deserialize('result/after_prediction.p')
# sample = 'SRR2603944_secondAligned.sortedByCoord.out.bed'
# jcmq.show_neoantigen_as_fasta(outdir='result/fasta',name='neoantigen_{}.fasta'.format(sample),stage=2,verbosity=1,contain_uid=True)
# snaf.remove_redundant('result/fasta/neoantigen_SRR2603944_secondAligned.sortedByCoord.out.bed.fasta','result/fasta/neoantigen_SRR2603944_secondAligned.sortedByCoord.out.bed_unique.fasta')
# snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot_9_10_mers_unique.fasta',
#                                   fa2_path='result/fasta/neoantigen_SRR2603944_secondAligned.sortedByCoord.out.bed_unique.fasta',outdir='result/fasta',
#                                   write_unique2=True,prefix='neoantigen_SRR2603944_secondAligned.sortedByCoord.out.bed_unique_')



# # get neoantigens whose junction are present in this cell line
# tcga = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/T_candidates/T_antigen_candidates_all.txt',sep='\t',index_col=0)
# hla = ['HLA-A*03:01','HLA-A*01:01','HLA-B*07:02','HLA-C*07:02']
# tcga = tcga.loc[tcga['hla'].isin(hla),:]
# df = snaf.remove_trailing_coord('altanalyze_output/ExpressionInput/counts.original.txt',sep='\t')
# # df = df.loc[df.mean(axis=1) > 20,:]
# uid_cl = set(df.index)
# mapping = df.iloc[:,0].to_dict()
# tcga_subset = tcga.loc[tcga['uid'].isin(uid_cl),:]
# tcga_subset['cell_line_count'] = tcga_subset['uid'].map(mapping).fillna(0).values

# # write and generate fasta
# tcga_subset.to_csv('tcga_melanoma_present.txt',sep='\t');sys.exit('stop')
# with open('search_db.fasta','w') as f:
#     for aa,sub_df in tcga_subset.groupby(by='peptide'):
#         f.write('>' + ','.join(list(set(sub_df['uid'].tolist()))) + '\n')
#         f.write('{}\n'.format(aa))

# polish fasta
snaf.remove_redundant('search_db.fasta','search_db_unique.fasta')
snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot_9_10_mers_unique.fasta',
                                  fa2_path='search_db_unique.fasta',outdir='.',
                                  write_unique2=True,prefix='search_db_unique_')  
                        

# configure maxquant
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/melanoma_ms/3784/search_db_unique_unique2.fasta']
inputs = [os.path.join('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/melanoma_ms/3784/raw','NCI_3784Mel_Epitope_Rep{}.raw'.format(i+1)) for i in range(5)]
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/melanoma_ms/3784/raw'
snaf.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)











