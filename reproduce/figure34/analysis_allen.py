#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
sys.path.insert(0,'/data/salomonis2/software')
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import re

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


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
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
db_dir = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'

snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='netMHCpan',software_path=netMHCpan_path)
surface.initialize(db_dir=db_dir)

# T antigen
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=50)
# hlas = [hla_string.split(',') for hla_string in df.columns.map(full_dict)]
# jcmq.run(hlas=hlas,outdir='./result')
# snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')

# 1. survival analysis
survival = pd.read_csv('all_patients_hla.txt',sep='\t',index_col=0)  
burden = pd.read_csv('result/burden_stage0.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1] 
survival = survival.loc[survival.index.isin(patient_dict.keys()),:]
reverse_patient_dict = {v:k for k,v in patient_dict.items()}
burden.index = burden.index.map(reverse_patient_dict).values
burden,burden_encode,burden_vc = snaf.survival_analysis(burden,survival,n=2,stratification_plot='result/stage0_stratify.pdf',survival_plot='result/stage0_survival.pdf',survival_duration='overall_survival',survival_event='dead')
sys.exit('stop')

# 2. embed neoantigen to umap
# snaf.downstream.analyze_neoantigens(freq_path='result/frequency_stage2_verbosity1_uid.txt',junction_path='result/burden_stage0.txt',total_samples=39,outdir='result',mers=None,fasta=False)
# snaf.run_dash_T_antigen(input_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/immunotherapy/allen/result/shared_vs_unique_neoantigen_all.txt',
#                         output_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/immunotherapy/allen/result')
# snaf.downstream.plot_umap_neoantigen(df_path='result/mer10_umap_embed_df.txt',outdir='result')


# survival correlation
# freq = 0.05
# tcga = pd.read_csv('survival_correlation/frequency_stage2_verbosity1_uid_TCGA_SKCM.txt',sep='\t',index_col=0)
# tcga.index = [item.replace(',','.').replace(':','.').replace('-','.') for item in tcga.index]
# tcga = tcga.loc[tcga['n_sample']>freq*472,:]
# tcga_dict = (tcga['n_sample']/472).to_dict()
# allen = pd.read_csv('survival_correlation/frequency_stage2_verbosity1_uid_van_allen.txt',sep='\t',index_col=0)
# allen.index = [item.replace(',','.').replace(':','.').replace('-','.') for item in allen.index]
# allen = allen.loc[allen['n_sample']>freq*39,:]
# allen_dict = (allen['n_sample']/39).to_dict()

# surv_corr_tcga = pd.read_csv('survival_correlation/SurvResult_dPSI_stage2_TCGA_SKCM.txt',sep='\t',index_col=0)
# surv_corr_tcga = surv_corr_tcga.loc[surv_corr_tcga.index.isin(tcga.index),:]
# surv_corr_allen = pd.read_csv('survival_correlation/SurvResult_dPSI_stage2_VanAllen.txt',sep='\t',index_col=0)
# surv_corr_allen = surv_corr_allen.loc[surv_corr_allen.index.isin(allen.index),:]

# surv_corr_merge = surv_corr_tcga.join(other=surv_corr_allen,how='inner',lsuffix='_tcga',rsuffix='_allen')
# surv_corr_merge = surv_corr_merge.loc[(surv_corr_merge['LRT_tcga']<0.05)&(surv_corr_merge['LRT_allen']<0.05),:]
# surv_corr_merge['freq_tcga'] = surv_corr_merge.index.map(tcga_dict).values
# surv_corr_merge['freq_allen'] = surv_corr_merge.index.map(allen_dict).values
# surv_corr_merge.to_csv('survival_correlation/freq_{}_result.txt'.format(freq),sep='\t')

df = pd.read_csv('survival_correlation/freq_0.05_result.txt',sep='\t',index_col=0)
fig,ax = plt.subplots()
sp = ax.scatter(df['Zscore_tcga'].values,df['Zscore_allen'].values,s=20,c=df['freq_tcga'].values)
for direc in ['top','right']:
    ax.spines[direc].set_visible(False)
for direc in ['bottom','left']:
    ax.spines[direc].set_position('zero')
ax.set_xlim((-6.01,6.01))
ax.set_ylim((-5.01,5.01))
ax.set_xticks(np.arange(-6,7,1))
ax.set_xticklabels(np.arange(-6,7,1))
ax.set_yticks(np.arange(-5,6,1))
ax.set_yticklabels(np.arange(-5,6,1))
ax.grid(alpha=0.5)
plt.colorbar(mappable=sp)
plt.savefig('survival_correlation/freq_0.05.pdf',bbox_inches='tight')
plt.close()



























