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
from scipy.stats import mannwhitneyu
from tqdm import tqdm
import statsmodels.stats.multitest as ssm



# # preprocess the dataframe
# df = pd.read_csv('./counts.TCGA-SKCM.txt',sep='\t',index_col=0)
# df.index = [item.split('=')[0] for item in df.index]
# df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]

# # filter to EventAnnotation file
# ee = [':'.join(item.split('|')[0].split(':')[1:]) for item in pd.read_csv('Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t')['UID']]
# df = df.loc[df.index.isin(set(ee)),:]


# # debug
# sample_to_hla = pd.read_csv('../HLA_genotype/sample_hla.txt',sep='\t',index_col=0)['hla'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# from functools import reduce
# combined_unique_hlas = reduce(lambda a,b:list(set(a+b)),hlas)
# combined_unique_hlas = snaf.snaf.hla_formatting(combined_unique_hlas,'netMHCpan_output','deepimmuno')
# df = snaf.binding.run_MHCflurry(['VREGVQISSR','LGARRSADQQ','WGSVREGVQI'],combined_unique_hlas)
# df = snaf.binding.run_MHCflurry(['VREGVQISSR','LGARRSADQQ','WGSVREGVQI'],combined_unique_hlas)


# # debug
# df = df.sample(n=100,axis=0)

# run SNAF
exon_table = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_Ensembl_exon_add_col.txt'
fasta = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_gene-seq-2000_flank.fa'
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
gtex_db = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/GTEx_junction_counts.h5ad'

snaf.initialize(exon_table=exon_table,fasta=fasta,software_path=netMHCpan_path,gtex_db=gtex_db,binding_method='netMHCpan')
# jcmq_skcm = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=50)
# jcmq_skcm.parallelize_run(kind=1)
# sample_to_hla = pd.read_csv('../HLA_genotype/sample_hla.txt',sep='\t',index_col=0)['hla'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq_skcm.parallelize_run(kind=3,hlas=hlas)
# jcmq_skcm.serialize(outdir='.',name='after_prediction.p')
# jcmq_skcm = snaf.JunctionCountMatrixQuery.deserialize(name='./after_prediction.p')
# jcmq_skcm.show_neoantigen_burden(outdir='.',name='burden_stage1.txt',stage=1,verbosity=1)
# jcmq_skcm.show_neoantigen_frequency(outdir='.',name='frequency_stage1.txt',stage=1,verbosity=1,plot=True,plot_name='frequency_stage1.pdf')


# # post-analysis
# jcmq = snaf.JunctionCountMatrixQuery.deserialize('after_prediction.p')
# burden = pd.read_csv('burden_stage3.txt',sep='\t',index_col=0)
# burden = burden.loc[jcmq.valid,:]
# sns.histplot(burden['mean'].values,binwidth=0.5)
# plt.savefig('junction_potential.pdf',bbox_inches='tight')
# plt.close()


'''patient analysis'''
# # 1. survival analysis
# survival = pd.read_csv('patient_analysis/TCGA-SKCM.survival.tsv',sep='\t',index_col=0)  # 463
# burden = pd.read_csv('burden_stage2.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
# burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
# burden = burden.loc[burden.index.isin(survival.index)]   # 458
# burden_encode = pd.Series(data=pd.qcut(burden,4,labels=['low','medium_low','medium_high','high']).to_list(),index=burden.index)
# be_vc = burden_encode.value_counts()
# sns.boxplot(data=burden.values);plt.savefig('patient_analysis/stratification.pdf',bbox_inches='tight');plt.close()
# high_group = burden_encode.loc[burden_encode=='high'].index.tolist()
# low_group = burden_encode.loc[burden_encode=='low'].index.tolist()
# high_os = survival.loc[high_group,['OS','OS.time']]
# low_os = survival.loc[low_group,['OS','OS.time']]
# fig,ax = plt.subplots()
# ax.set_ylim(0,1)
# for df in [low_os,high_os]:
#     kmf = KaplanMeierFitter()
#     kmf.fit(df['OS.time'],df['OS'])
#     kmf.plot_survival_function(ax=ax,ci_show=False,at_risk_counts=False)
# current_handles,current_labels = ax.get_legend_handles_labels()
# new_labels = ['low_burden','high_burden']
# ax.legend(current_handles,new_labels)
# results = logrank_test(low_os['OS.time'],high_os['OS.time'],low_os['OS'],high_os['OS'])
# ax.text(x=1000,y=0.05,s='Log-rank test: p-value is {:.2f}'.format(results.p_value),weight='bold')
# plt.savefig('patient_analysis/survial_high_low.pdf',bbox_inches='tight');plt.close()

# 2. mutation analysis
# mutation = pd.read_csv('patient_analysis/TCGA-SKCM.mutect2_snv.tsv',sep='\t',index_col=0)  # 467 samples have mutations
# mutation = mutation.loc[mutation['filter']=='PASS',:]
# burden = pd.read_csv('burden_stage2.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
# burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
# burden = burden.loc[burden.index.isin(mutation.index.unique())]   # 465
# with open('patient_analysis/mutation.txt','w') as f:
#     f.write('mutation_gene\tn_samples\tpval\n')
#     for gene in tqdm(mutation['gene'].unique()):
#         yes_samples = mutation.loc[mutation['gene'] == gene,:].index.unique().tolist()
#         if len(yes_samples) < 10:
#             continue
#         burden_df = burden.to_frame()
#         burden_df.columns = ['burden']
#         burden_df['mutation_{}'.format(gene)] = [True if sample in set(yes_samples) else False for sample in burden_df.index]
#         x = burden_df.loc[burden_df['mutation_{}'.format(gene)],'burden'].values
#         y = burden_df.loc[~(burden_df['mutation_{}'.format(gene)]),'burden'].values
#         u,p = mannwhitneyu(x=x,y=y)
#         f.write('{}\t{}\t{}\n'.format(gene,len(yes_samples),p))
# asso = pd.read_csv('patient_analysis/mutation.txt',sep='\t',index_col=0)
# results = ssm.multipletests(asso['pval'].values,alpha=0.2,method='fdr_bh')
# asso['adjp'] = results[1]
# asso.to_csv('patient_analysis/mutation_adjp.txt',sep='\t')
# for gene in ['KCNJ5']:
#     yes_samples = mutation.loc[mutation['gene'] == gene,:].index.unique().tolist()
#     burden_df = burden.to_frame()
#     burden_df.columns = ['burden']
#     burden_df['mutation_{}'.format(gene)] = [True if sample in set(yes_samples) else False for sample in burden_df.index]
#     x = burden_df.loc[burden_df['mutation_{}'.format(gene)],'burden'].values
#     y = burden_df.loc[~(burden_df['mutation_{}'.format(gene)]),'burden'].values
#     u,p = mannwhitneyu(x=x,y=y)
#     fig,ax = plt.subplots()
#     sns.boxplot(data=burden_df,x='mutation_{}'.format(gene),y='burden',ax=ax,width=0.5)
#     ax.text(x=0.5,y=0.5,s='mannwhitney p={}'.format(p),weight='bold')
#     plt.savefig('patient_analysis/mutation_plot/{}_diff.pdf'.format(gene),bbox_inches='tight');plt.close()


'''junction analysis'''
# 1. tumor specificity
# jcmq_skcm = snaf.JunctionCountMatrixQuery.deserialize(name='./after_prediction.p')
burden = pd.read_csv('burden_stage3.txt',sep='\t',index_col=0)
nonzero = np.count_nonzero(burden.values,axis=1)/burden.shape[1]
burden = burden.loc[nonzero>0,:]
burden = burden.iloc[:-1,:]
# bayes
# bayes_col = []
# for uid in tqdm(burden.index):
#     sigma = snaf.gtex.accurate_tumor_specificity(uid,'bayesian')
#     bayes_col.append(sigma)
# burden['bayes_sigma'] = sigma_col
# crude mean
mean_col = []
for uid in tqdm(burden.index):
    _,sigma = snaf.gtex.crude_tumor_specificity(uid,50)
    mean_col.append(sigma)
burden['mean_gtex'] = mean_col
# mle
sigma_col = []
for uid in tqdm(burden.index):
    sigma = snaf.gtex.accurate_tumor_specificity(uid,'mle')
    sigma_col.append(sigma)
burden['mle_sigma'] = sigma_col
burden.loc[:,['mean','mean_gtex','mle_sigma']].to_csv('junction_analysis/specificity.txt',sep='\t')





