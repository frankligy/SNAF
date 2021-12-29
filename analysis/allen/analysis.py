#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
import snaf
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy.stats import mannwhitneyu,pearsonr,spearmanr
from tqdm import tqdm
import statsmodels.stats.multitest as ssm
import re


# get samples_hla file
rna = pd.read_csv('srr_list.txt',sep='\t',index_col=0,header=None)
pat = re.compile(r'Pat\d{2,3}')
patient = []
for item in rna.index:
    if not re.search(pat,item):
        raise Exception('check,stop')
    else:
        patient.append(re.search(pat,item).group(0))
rna['patient'] = patient
rna.set_index(keys='patient',inplace=True)
rna['full_name'] = [item + '_1_' + item + '_2.Aligned.sortedByCoord.out.bed' for item in rna[1]]
patient_dict = rna['full_name'].to_dict()   # {Pat40:SRR_1_SRR_2.Aligned..s}

hla = pd.read_csv('all_patients_hla.txt',sep='\t',index_col=0).loc[:,['hla.a1','hla.a2','hla.b1','hla.b2','hla.c1','hla.c2']]
hla = hla.loc[hla.index.isin(patient_dict.keys()),:]
hla_dict = {}      # {Pat40:hla_string}
for row in hla.iterrows():
    hla_dict[row[0]] = ','.join(snaf.snaf.hla_formatting(row[1].tolist(),pre_type='netMHCpan_input',post_type='netMHCpan_output'))

full_dict = {}   # {SRR_1...:hla_string}
for p,f in patient_dict.items():
    full_dict[f] = hla_dict[p]



# preprocess the dataframe
df = pd.read_csv('./counts.original.txt',sep='\t',index_col=0)
df.index = [item.split('=')[0] for item in df.index]
df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]  # 668986 rows x 39 columns


# filter to EventAnnotation file
ee = [':'.join(item.split('|')[0].split(':')[1:]) for item in pd.read_csv('Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t')['UID']]
df = df.loc[df.index.isin(set(ee)),:]   # 49114 rows x 39 columns


# run SNAF
exon_table = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_Ensembl_exon_add_col.txt'
fasta = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_gene-seq-2000_flank.fa'
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
gtex_db = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/GTEx_junction_counts.h5ad'

# snaf.initialize(exon_table=exon_table,fasta=fasta,software_path=netMHCpan_path,gtex_db=gtex_db,binding_method='netMHCpan')
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=50)
# jcmq.parallelize_run(kind=1)
# hlas = [hla_string.split(',') for hla_string in df.columns.map(full_dict)]
# jcmq.parallelize_run(kind=3,hlas=hlas)
# jcmq.serialize(outdir='.',name='after_prediction.p')
# jcmq = snaf.JunctionCountMatrixQuery.deserialize(name='./after_prediction.p')
# burden_stage0 = jcmq.cond_df.astype('int8')
# burden_stage0.loc['burden'] = burden_stage0.sum(axis=0).values
# burden_stage0['mean'] = burden_stage0.mean(axis=1).values
# burden_stage0.to_csv('burden_stage0.txt',sep='\t')
# for stage in [3,2,1]:
#     jcmq.show_neoantigen_burden(outdir='.',name='burden_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False)
#     jcmq.show_neoantigen_frequency(outdir='.',name='frequency_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False,plot=True,plot_name='frequency_stage{}.pdf'.format(stage))
# jcmq.show_neoantigen_frequency(outdir='.',name='frequency_stage3_verbosity1_uid.txt',stage=3,verbosity=1,contain_uid=True,plot=False)


'''patient analysis'''
# 1. survival analysis
survival = pd.read_csv('all_patients_hla.txt',sep='\t',index_col=0)  
burden = pd.read_csv('burden_stage2.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1] 
survival = survival.loc[survival.index.isin(patient_dict.keys()),:]
reverse_patient_dict = {v:k for k,v in patient_dict.items()}
burden.index = burden.index.map(reverse_patient_dict).values
# burden,burden_encode,burden_vc = snaf.survival_analysis(burden,survival,n=2,stratification_plot='patient_analysis/stage2_stratify.pdf',survival_plot='patient_analysis/stage2_survival.pdf',survival_duration='overall_survival',survival_event='dead')
# burden,burden_encode,burden_vc = snaf.survival_analysis(burden,survival,n=2,stratification_plot='patient_analysis/stage0_stratify.pdf',survival_plot='patient_analysis/stage0_progression.pdf',survival_duration='progression_free',survival_event='progression')


# 2. all burden
burden0 = pd.read_csv('burden_stage0.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]
burden1 = pd.read_csv('burden_stage1.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]
burden2 = pd.read_csv('burden_stage2.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]
burden3 = pd.read_csv('burden_stage3.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]
burden = pd.concat([burden0,burden1,burden2,burden3],axis=1)
burden.columns = ['burden0','burden1','burden2','burden3']
burden.to_csv('patient_analysis/burden_dynamics.txt',sep='\t')
sys.exit('stop')

'''neoantigen analysis'''
# 1. physicalchemical properties relate to occurence?
# freq3 = pd.read_csv('neoantigen_analysis/frequency_stage3_verbosity1_uid.txt',sep='\t',index_col=0)
# burden0 = pd.read_csv('burden_stage0.txt',sep='\t',index_col=0)['mean'].to_dict()
# freq3['uid'] = [item.split(',')[-1] for item in freq3.index]
# freq3.index = [item.split(',')[0] for item in freq3.index]
# freq3['burden0'] = freq3['uid'].map(burden0).values
# freq3['expected'] = freq3['burden0'] * 39
# freq3.drop(columns='samples',inplace=True)
# identity = []
# for o,e in zip(freq3['n_sample'].values,freq3['expected'].values):
#     if o < 0.1 * e:
#         identity.append('low')
#     elif o > 0.9 * e:
#         identity.append('high')
#     else:
#         identity.append('medium')
# freq3['identity'] = identity
# freq3 = freq3.loc[np.logical_not(freq3.index.duplicated()),:]
# freq3 = freq3.loc[freq3['identity']!='medium',:]
# freq3 = freq3.loc[freq3['burden0']>0.1,:]
# freq3.to_csv('neoantigen_analysis/df_test_app.txt',sep='\t')
# sns.regplot(freq3['burden0'].values,freq3['n_sample'])
# plt.savefig('neoantigen_analysis/freq3_correlation_filter.pdf',bbox_inches='tight');plt.close()

snaf.run_dash_app(intpath='/data/salomonis2/LabFiles/Frank-Li/neoantigen/immunotherapy/allen/neoantigen_analysis/df_test_app.txt',
                  remove_cols=['uid'],host='bmi-460g9-20.chmcres.cchmc.org')

























