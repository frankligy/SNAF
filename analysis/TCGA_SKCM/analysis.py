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
    #jcmq_skcm.show_neoantigen_burden(outdir='.',name='burden_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False)
    # jcmq_skcm.show_neoantigen_frequency(outdir='.',name='frequency_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False,plot=True,plot_name='frequency_stage{}.pdf'.format(stage))
# jcmq_skcm.show_neoantigen_frequency(outdir='.',name='frequency_stage3_verbosity1_uid.txt',stage=3,verbosity=1,contain_uid=True,plot=False)



'''patient analysis'''
# 1. survival analysis
# survival = pd.read_csv('patient_analysis/TCGA-SKCM.survival.tsv',sep='\t',index_col=0)  # 463
# burden = pd.read_csv('burden_stage0.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
# burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
# burden,burden_encode,burden_vc = snaf.survival_analysis(burden,survival,stratification_plot='patient_analysis/stage0_stratify.pdf',survival_plot='patient_analysis/stage0_survival.pdf')

# 2. mutation analysis
# mutation = pd.read_csv('patient_analysis/TCGA-SKCM.mutect2_snv.tsv',sep='\t',index_col=0)  # 467 samples have mutations
# mutation = mutation.loc[mutation['filter']=='PASS',:]
# burden = pd.read_csv('burden_stage0.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
# burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
# # snaf.mutation_analysis(mode='compute',burden=burden,mutation=mutation,output='patient_analysis/stage3_mutation.txt')
# snaf.mutation_analysis(mode='plot',burden=burden,mutation=mutation,output='patient_analysis/stage3_mutation_CAMKK2.pdf',genes_to_plot=['CAMKK2'])

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


'''neoantigen analysis'''
# 1. physicalchemical properties relate to occurence?
# freq3 = pd.read_csv('frequency_stage3_verbosity1_uid.txt',sep='\t',index_col=0)
# burden0 = pd.read_csv('burden_stage0.txt',sep='\t',index_col=0)['mean'].to_dict()
# freq3['uid'] = [item.split(',')[-1] for item in freq3.index]
# freq3.index = [item.split(',')[0] for item in freq3.index]
# freq3['burden0'] = freq3['uid'].map(burden0).values
# freq3['expected'] = freq3['burden0'] * 472
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

snaf.run_dash_app(intpath='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/neoantigen_analysis/df_test_app.txt',
                  remove_cols=['uid'],host='bmi-460g9-10.chmcres.cchmc.org')

























