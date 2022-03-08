#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm


# preprocess the dataframe
df = pd.read_csv('./counts.TCGA-SKCM.txt',sep='\t',index_col=0)
df.index = [item.split('=')[0] for item in df.index]
df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]

# filter to EventAnnotation file
ee = [':'.join(item.split('|')[0].split(':')[1:]) for item in pd.read_csv('Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t')['UID']]
df = df.loc[df.index.isin(set(ee)),:]


# run SNAF
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
db_dir = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'

snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='netMHCpan',software_path=netMHCpan_path)
surface.initialize(db_dir=db_dir)

# # gtex check
# snaf.gtex_visual_combine(query='ENSG00000101003:E7.2-E9.2',norm=True,tumor=df)


'''B cell neoantigen'''
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=30)
# neojunctions = jcmq.valid
# membrane_uid = surface.filter_to_membrane_protein(neojunctions)
# membrane_uid = [(uid,snaf.gtex.accurate_tumor_specificity(uid,method='mean'),jcmq.get_neojunction_info(uid)[0],jcmq.get_neojunction_info(uid)[1]) for uid in membrane_uid]
# surface.single_run(membrane_uid,2,True,'/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',serialize=True)
# c_candidates,c_further = surface.process_results('single_run_surface_antigen.p',strigency=3)
# print(c_candidates,c_further)

# sa = snaf.surface.individual_check('ENSG00000176204:E5.2-E6.1',tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',indices=[None],fragments=[None])
# surface.run_dash_prioritizer(pkl='single_run_surface_antigen.p',candidates='candidates.txt',python_executable='/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7')
# snaf.gtex_viewer.gtex_visual(query='ENSG00000090339:E4.3-E4.5',norm=True)
# snaf.gtex_viewer.gtex_visual(query='ENSG00000090339:E4.3-E4.5',norm=False)


'''T cell neoantigen'''
#### Step1: Running the program
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=50)
# sample_to_hla = pd.read_csv('../HLA_genotype/sample_hla.txt',sep='\t',index_col=0)['hla'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq.run(hlas=hlas,outdir='./result')
# snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')


### Step2: necessary secondary results
# df = pd.read_csv('./result/frequency_stage3_verbosity1_uid.txt',sep='\t',index_col=0)
# snaf.downstream.add_gene_symbol_frequency_table(df=df).to_csv('./result/frequency_stage3_verbosity1_uid_gene_symbol.txt',sep='\t')
# jcmq = snaf.JunctionCountMatrixQuery.deserialize(name='./result/after_prediction.p')
# jcmq.visualize(uid='ENSG00000198034:E8.4-E9.1',sample='TCGA-WE-A8ZT-06A-11R-A37K-07.bed',outdir='./result')


### Step3: downstream analysis (patient level and neoantigen level)

'''patient analysis'''
# 1. survival analysis
# survival = pd.read_csv('patient_analysis/TCGA-SKCM.survival.tsv',sep='\t',index_col=0)  # 463
# burden = pd.read_csv('burden_stage2.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
# burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
# burden,burden_encode,burden_vc = snaf.survival_analysis(burden,survival,n=2,stratification_plot='patient_analysis/stage2_stratify.pdf',survival_plot='patient_analysis/stage2_survival.pdf')
# burden_encode.to_csv('to_anu/stage2_stratify.txt',sep='\t',header=None)

# 2. mutation analysis
# mutation = pd.read_csv('patient_analysis/TCGA-SKCM.mutect2_snv.tsv',sep='\t',index_col=0)  # 467 samples have mutations
# mutation = mutation.loc[mutation['filter']=='PASS',:]
# burden = pd.read_csv('burden_stage3.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
# burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
# # snaf.mutation_analysis(mode='compute',burden=burden,mutation=mutation,output='patient_analysis/stage3_mutation.txt')
# snaf.mutation_analysis(mode='plot',burden=burden,mutation=mutation,output='patient_analysis/stage3_mutation_CAMKK2.pdf',genes_to_plot=['CAMKK2'])


'''neoantigen analysis'''
# 1. physicalchemical properties relate to occurence?
# freq = pd.read_csv('to_anu/frequency_stage2_verbosity1_uid.txt',sep='\t',index_col=0)
# burden0 = pd.read_csv('burden_stage0.txt',sep='\t',index_col=0)['mean'].to_dict()
# freq['uid'] = [item.split(',')[-1] for item in freq.index]
# freq.index = [item.split(',')[0] for item in freq.index]
# freq['burden0'] = freq['uid'].map(burden0).values
# freq['expected'] = freq['burden0'] * 472
# freq.drop(columns='samples',inplace=True)
# identity = []
# for o,e in zip(freq['n_sample'].values,freq['expected'].values):
#     if o < 0.1 * e:
#         identity.append('low')
#     elif o > 0.9 * e:
#         identity.append('high')
#     else:
#         identity.append('medium')
# freq['identity'] = identity
# freq = freq.loc[np.logical_not(freq.index.duplicated()),:]
# freq = freq.loc[freq['identity']!='medium',:]
# freq = freq.loc[freq['burden0']>0.1,:]
# freq['length'] = [len(item) for item in freq.index]
# freq_9mer = freq.loc[freq['length']==9,:]
# freq_10mer = freq.loc[freq['length']==10,:]
# freq.to_csv('to_anu/neoantigen_common_unique.txt',sep='\t')
# freq_9mer.to_csv('neoantigen_analysis/df_test_app_mer9.txt',sep='\t')
# freq_10mer.to_csv('neoantigen_analysis/df_test_app_mer10.txt',sep='\t')
# sns.regplot(freq['burden0'].values,freq['n_sample'])
# plt.savefig('neoantigen_analysis/freq_correlation_filter.pdf',bbox_inches='tight');plt.close()

# snaf.run_dash_app(intpath='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/neoantigen_analysis/df_test_app.txt',
#                   remove_cols=['uid'],host='bmi-460g9-10.chmcres.cchmc.org')

# # discriminative motif analysis
# df = pd.read_csv('neoantigen_analysis/df_test_app_mer9.txt',sep='\t',index_col=0)
# with open('neoantigen_analysis/mer9_high.fasta','w') as f1, open('neoantigen_analysis/mer9_low.fasta','w') as f2:
#     for identity,sub_df in df.groupby(by='identity'):
#         if identity == 'high':
#             for item in sub_df.index:
#                 f1.write('>{}\n{}\n'.format(item,item))
#         else:
#             for item in sub_df.index:
#                 f2.write('>{}\n{}\n'.format(item,item))          




























