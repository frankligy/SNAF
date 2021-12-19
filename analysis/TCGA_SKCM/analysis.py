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
from ast import literal_eval



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

# snaf.initialize(exon_table=exon_table,fasta=fasta,software_path=netMHCpan_path,gtex_db=gtex_db,binding_method='netMHCpan')
# jcmq_skcm = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=50)
# jcmq_skcm.parallelize_run(kind=1)
# sample_to_hla = pd.read_csv('../HLA_genotype/sample_hla.txt',sep='\t',index_col=0)['hla'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq_skcm.parallelize_run(kind=3,hlas=hlas)
# jcmq_skcm.serialize(outdir='.',name='after_prediction.p')
# jcmq_skcm = snaf.JunctionCountMatrixQuery.deserialize(name='./after_prediction.p')
# burden_stage0 = jcmq_skcm.cond_df.astype('int8')
# burden_stage0.loc['burden'] = burden_stage0.sum(axis=0).values
# burden_stage0['mean'] = burden_stage0.mean(axis=1).values
# burden_stage0.to_csv('burden_stage0.txt',sep='\t')
# for stage in [3,2,1]:
#     jcmq_skcm.show_neoantigen_burden(outdir='.',name='burden_stage{}.txt'.format(stage),stage=stage,verbosity=1)
#     jcmq_skcm.show_neoantigen_frequency(outdir='.',name='frequency_stage{}.txt'.format(stage),stage=stage,verbosity=1,plot=True,plot_name='frequency_stage{}.pdf'.format(stage))
# for stage in [3,2,1]:
#     jcmq_skcm.show_neoantigen_frequency(outdir='.',name='frequency_stage{}_uid.txt'.format(stage),stage=stage,verbosity=1,contain_uid=True,plot=False)

# # post-analysis
# jcmq = snaf.JunctionCountMatrixQuery.deserialize('after_prediction.p')
# burden = pd.read_csv('burden_stage3.txt',sep='\t',index_col=0)
# burden = burden.loc[jcmq.valid,:]
# sns.histplot(burden['mean'].values,binwidth=0.5)
# plt.savefig('junction_potential.pdf',bbox_inches='tight')
# plt.close()


'''patient analysis'''
# 1. survival analysis
# survival = pd.read_csv('patient_analysis/TCGA-SKCM.survival.tsv',sep='\t',index_col=0)  # 463
# burden = pd.read_csv('burden_stage0.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
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

# # 2. mutation analysis
# mutation = pd.read_csv('patient_analysis/TCGA-SKCM.mutect2_snv.tsv',sep='\t',index_col=0)  # 467 samples have mutations
# mutation = mutation.loc[mutation['filter']=='PASS',:]
# burden = pd.read_csv('burden_stage0.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
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
# results = ssm.multipletests(asso['pval'].values,alpha=0.1,method='fdr_bh')
# asso['adjp'] = results[1]
# asso.to_csv('patient_analysis/mutation_adjp.txt',sep='\t')
# for gene in ['CAMKK2']:
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

# 3. differential gene analysis
# count = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/patient_analysis/de/TCGA-SKCM.htseq_counts.tsv',sep='\t',index_col=0)  # (60488, 472)
# count.index = [ensg.split('.')[0] for ensg in count.index]   # the decimal part indicate how many times this model has been changed in Ensembl
# count = count.loc[np.logical_not(count.index.duplicated()),:] # (60488, 472)

# def ensemblgene_to_symbol(query,species):
#     '''
#     Examples::
#         from sctriangulate.preprocessing import GeneConvert
#         converted_list = GeneConvert.ensemblgene_to_symbol(['ENSG00000010404','ENSG00000010505'],species='human')
#     '''
#     # assume query is a list, will also return a list
#     import mygene
#     mg = mygene.MyGeneInfo()
#     out = mg.querymany(query,scopes='ensemblgene',fileds='symbol',species=species,returnall=True,as_dataframe=True,df_index=True)
#     result = out['out']['symbol'].fillna('unknown_gene').tolist()
#     try:
#         assert len(query) == len(result)
#     except AssertionError:    # have duplicate results
#         print('having duplicated results, only take unique ones')
#         df = out['out']
#         df_unique = df.loc[np.logical_not(df.index.duplicated()),:]
#         df_unique = df_unique.loc[query,:]
#         result = df_unique['symbol'].fillna('unknown_gene').tolist()
#     return result

# gene_symbol = ensemblgene_to_symbol(query=count.index.tolist(),species='human')
# count.index = gene_symbol
# count = count.loc[count.index!='unknown_gene',:]  # (40443, 472)
# count = count.loc[np.logical_not(count.index.duplicated()),:]  # (39271, 472)
# count = pd.DataFrame(data=np.clip(np.round_(np.exp2(count.values)),a_min=1,a_max=None) - 1,index=count.index,columns=count.columns)
# count.to_csv('patient_analysis/de/gene_symbol_count_matrix.txt',sep='\t')

# count = pd.read_csv('patient_analysis/de/gene_symbol_count_matrix.txt',sep='\t',index_col=0)
# survival = pd.read_csv('patient_analysis/TCGA-SKCM.survival.tsv',sep='\t',index_col=0)  # 463
# burden = pd.read_csv('burden_stage3.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
# burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
# burden = burden.loc[burden.index.isin(survival.index)]   # 458
# burden_encode = pd.Series(data=pd.qcut(burden,4,labels=['low','medium_low','medium_high','high']).to_list(),index=burden.index)
# burden_encode = burden_encode.loc[(burden_encode=='low')|(burden_encode=='high')].to_frame(name='group')
# count_s = count.loc[:,burden_encode.index]
# count_s.to_csv('patient_analysis/de/deseq2_count.txt',sep='\t')
# burden_encode.to_csv('patient_analysis/de/deseq2_meta.txt',sep='\t')

# deseq_de = pd.read_csv('patient_analysis/de/deseq2_results.txt',sep='\t',index_col=0)
# high_gene = deseq_de.loc[(deseq_de['log2FoldChange']>2)&(deseq_de['padj']<0.01),:].index.tolist()
# low_gene = deseq_de.loc[(deseq_de['log2FoldChange']<-2)&(deseq_de['padj']<0.01),:].index.tolist()


# 4. dynamic analysis
# burden0 = pd.read_csv('burden_stage0.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]
# burden1 = pd.read_csv('burden_stage1.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]
# burden2 = pd.read_csv('burden_stage2.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]
# burden3 = pd.read_csv('burden_stage3.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]
# burden = pd.concat([burden0,burden1,burden2,burden3],axis=1)
# burden.columns = ['burden0','burden1','burden2','burden3']
# burden.to_csv('patient_analysis/burden_dynamics.txt',sep='\t')











'''junction analysis'''
# 1. tumor specificity
# jcmq_skcm = snaf.JunctionCountMatrixQuery.deserialize(name='./after_prediction.p')
# burden = pd.read_csv('burden_stage3.txt',sep='\t',index_col=0)
# nonzero = np.count_nonzero(burden.values,axis=1)/burden.shape[1]
# burden = burden.loc[nonzero>0,:]
# burden = burden.iloc[:-1,:]
# bayes, takes so long, just sanity check
# for uid in tqdm(burden.index[:5]):
#     sigma = snaf.gtex.accurate_tumor_specificity(uid,'bayesian')
# # crude mean
# mean_col = []
# for uid in tqdm(burden.index):
#     _,sigma = snaf.gtex.crude_tumor_specificity(uid,50)
#     mean_col.append(sigma)
# burden['mean_gtex'] = mean_col
# # mle
# sigma_col = []
# for uid in tqdm(burden.index):
#     sigma = snaf.gtex.accurate_tumor_specificity(uid,'mle')
#     sigma_col.append(sigma)
# burden['mle_sigma'] = sigma_col
# burden.loc[:,['mean','mean_gtex','mle_sigma']].to_csv('junction_analysis/specificity.txt',sep='\t')
# specificity = pd.read_csv('junction_analysis/specificity.txt',sep='\t',index_col=0)
# r,p = pearsonr(specificity['mean'].values,specificity['mean_gtex'].values)  
# r,p = spearmanr(specificity['mean'].values,specificity['mean_gtex'].values)  
# sns.regplot(specificity['mean'].values,specificity['mean_gtex'].values,scatter_kws={'s':3},line_kws={'color':'k','linewidth':3});plt.savefig('junction_analysis/specificity_to_mean_gtex.pdf',bbox_inches='tight');plt.close()
# r,p = pearsonr(specificity['mean'].values,specificity['mle_sigma'].values)  
# r,p = spearmanr(specificity['mean'].values,specificity['mle_sigma'].values)  
# sns.regplot(specificity['mean'].values,specificity['mle_sigma'].values,scatter_kws={'s':3},line_kws={'color':'k','linewidth':3});plt.savefig('junction_analysis/specificity_to_mle_sigma.pdf',bbox_inches='tight');plt.close()
# r,p = pearsonr(specificity['mean_gtex'].values,specificity['mle_sigma'].values)  # r=-0.35, p=0
# r,p = spearmanr(specificity['mean_gtex'].values,specificity['mle_sigma'].values) # r=-0.34, p=0
# sns.regplot(specificity['mean_gtex'].values,specificity['mle_sigma'].values,scatter_kws={'s':3},line_kws={'color':'k','linewidth':3});plt.savefig('junction_analysis/specificity_mean_gtex_vs_mle_sigma.pdf',bbox_inches='tight');plt.close()

# 2. number of neoantigens dynamics at different stages
# burden0 = pd.read_csv('burden_stage0.txt',sep='\t',index_col=0).loc[:,['mean']]
# burden1 = pd.read_csv('burden_stage1.txt',sep='\t',index_col=0).loc[:,['mean']]
# burden2 = pd.read_csv('burden_stage2.txt',sep='\t',index_col=0).loc[:,['mean']]
# burden3 = pd.read_csv('burden_stage3.txt',sep='\t',index_col=0).loc[:,['mean']]
# burden = pd.concat([burden0,burden1,burden2,burden3],axis=1)
# burden.columns = ['mean0','mean1','mean2','mean3']
# burden.to_csv('junction_analysis/burden_dynamics.txt',sep='\t')



'''neoantigen analysis'''
# 1. physicalchemical properties relate to occurence?
freq3_uid = pd.read_csv('frequency_stage3_uid.txt',sep='\t',index_col=0)
burden0 = pd.read_csv('burden_stage0.txt',sep='\t',index_col=0)['mean'].to_dict()
freq3_uid['uid'] = [item.split(',')[1] for item in freq3_uid.index]
freq3_uid.index = [item.split(',')[0] for item in freq3_uid.index]
freq3_uid['burden0'] = freq3_uid['uid'].map(burden0).values
freq3_uid['expected'] = freq3_uid['burden0'] * 472
identity = []
for o,e in zip(freq3_uid['n_sample'].values,freq3_uid['expected'].values):
    if o < 0.1 * e:
        identity.append('low')
    elif o > 0.9 * e:
        identity.append('high')
    else:
        identity.append('medium')
freq3_uid['identity'] = identity
freq3_uid = freq3_uid.loc[np.logical_not(freq3_uid.index.duplicated()),:]
# sns.regplot(freq3_uid['burden0'].values,freq3_uid['n_sample'])
# plt.savefig('neoantigen_analysis/freq3_correlation_raw.pdf',bbox_inches='tight');plt.close()
freq3_uid = freq3_uid.loc[freq3_uid['identity']!='medium',:]
# sns.regplot(freq3_uid['burden0'].values,freq3_uid['n_sample'])
# plt.savefig('neoantigen_analysis/freq3_correlation_filter.pdf',bbox_inches='tight');plt.close()

# pearsonr(freq3_uid['n_sample'].values,freq3_uid['burden0'].values)  # 0.66,0
# spearmanr(freq3_uid['n_sample'].values,freq3_uid['burden0'].values) # 0.73,0
freq3_uid.to_csv('neoantigen_analysis/freq3_source.txt',sep='\t')


freq3 = freq3_uid.loc[freq3_uid['burden0']>0.1,:]

after_pca = np.loadtxt('after_pca.txt')
def aaindex(peptide,after_pca):
    amino = 'ARNDCQEGHILKMFPSTWYV-'
    matrix = np.transpose(after_pca)   # [12,21]
    encoded = np.empty([len(peptide), 12])  # (seq_len,12)
    for i in range(len(peptide)):
        query = peptide[i]
        if query == 'X': query = '-'
        query = query.upper()
        encoded[i, :] = matrix[:, amino.index(query)]
    return encoded.flatten()

freq3['length'] = [len(item) for item in freq3.index]
freq3_mer9 = freq3.loc[freq3['length']==9,:]
freq3_mer10 = freq3.loc[freq3['length']==10,:]

# 9mer
mer9_encoded = np.empty([freq3_mer9.shape[0],9*12])
for i,pep in enumerate(freq3_mer9.index):
    mer9_encoded[i,:] = aaindex(pep,after_pca)
from sklearn.decomposition import PCA
from umap import UMAP
model = PCA()
model.fit(mer9_encoded)
accumulated_variance = np.cumsum(model.explained_variance_ratio_)  # 54 PCs alraedy explain 90%
pca_scoring = model.transform(mer9_encoded)[:,:55]
reducer = UMAP(random_state=42,min_dist=0.5)
embedding = reducer.fit_transform(pca_scoring)

df = pd.DataFrame(data={'umap_x':embedding[:,0],'umap_y':embedding[:,1],
                        'n_sample':freq3_mer9['n_sample'].values,'burden0':freq3_mer9['burden0'].values,
                        'identity':freq3_mer9['identity'],'uid':freq3_mer9['uid']},index=freq3_mer9.index)              
# fig,ax = plt.subplots()
# df_high = df.loc[df['identity']=='high',:]
# df_low = df.loc[df['identity']=='low',:]
# ax.scatter(x=df_low['umap_x'].values,y=df_low['umap_y'].values,c='lightgrey',s=0.5)
# ax.scatter(x=df_high['umap_x'].values,y=df_high['umap_y'].values,c='red',s=0.5)
# plt.savefig('neoantigen_analysis/mer9_umap_high_low.pdf',bbox_inches='tight');plt.close()

# test plotly
import plotly.graph_objects as go
fig = go.Figure(data=go.Scatter(x=df['umap_x'],y=df['umap_y'],mode='markers',
                                marker_color=['red' if item == 'high' else 'lightgrey' for item in df['identity']],
                                text=df.index.values,customdata=df['uid'].values,
                                hovertemplate='%{text}<br>%{customdata}'))
fig.update_layout(title='test',margin=dict(l=300,r=300),plot_bgcolor='rgba(0,0,0,0)')
fig.write_html('neoantigen_analysis/mer9_umap_high_low_interact.html',include_plotlyjs='cdn')





















