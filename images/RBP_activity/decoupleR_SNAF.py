#!/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/decoupler_env/bin/python3.8

import decoupler as dc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os,sys
from tqdm import tqdm
from adjustText import adjust_text

import matplotlib as mpl

# for publication ready figure
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# # preprocessing to get necessary input files, code derived from https://github.com/frankligy/splice-inferelator/blob/main/scripts/prepare_inference.py
# def disentangle_uid(uid,add_stable_id=False,add_symbol=False):
#     fg,bg = uid.split('|')
#     gene_symbol = fg.split(':')[0]
#     stable_id = fg.split(':')[1]
#     exons = ':'.join(fg.split(':')[2:])
#     output = exons
#     if add_stable_id:
#         output = stable_id + ':' + output
#     if add_symbol:
#         output = gene_symbol + ':' + output
#     return output

# def median_impute(x):
#     med = np.ma.median(np.ma.masked_invalid(x.values))
#     result = np.nan_to_num(x.values,nan=med)
#     return result

# ee = pd.read_csv('../TCGA_melanoma/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t',index_col='UID').iloc[:,10:]
# ee.index = [disentangle_uid(item,add_symbol=True) for item in ee.index]
# prop = 1- np.count_nonzero(ee.isna().values,axis=1)/ee.shape[1]  # 101974 rows x 472 columns
# ee = ee.loc[prop>0.75,:]   # 59134 rows x 472 columns
# original_columns = ee.columns
# ee = ee.apply(func=median_impute,axis=1,result_type='expand')
# ee.columns = original_columns
# cv = ee.apply(func=lambda x:x.values.std()/x.values.mean(),axis=1,result_type='reduce').values
# ee = ee.loc[cv>0.1,:]
# ee.columns = original_columns   # 34683 rows x 472 columns
# ee.T.to_csv('expr.tsv',sep='\t')


# # test on TCGA melanoma
# # expr: sample * event
# # prior: event * sf
expr_path = 'expr.tsv'
prior_path = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/HepG2_inference/prior.tsv'
root_path = '.'
expr = pd.read_csv(expr_path,sep='\t',index_col=0)  # 472 rows x 34683 columns
prior = pd.read_csv(prior_path,sep='\t',index_col=0)  # 73157 rows x 221 columns
common_events = list(set(expr.columns).intersection(set(prior.index)))
expr = expr.loc[:,common_events]   # 472 rows x 17836 columns
prior = prior.loc[common_events,:]  # 17836 rows x 221 columns
net = prior.stack().reset_index()
net.columns = ['target','source','weight']
net = net.loc[net['weight']!=0,:]

# # explore expr a bit
# sns.histplot(np.count_nonzero(expr.values,axis=1),bins=100)
# plt.savefig(os.path.join(root_path,'expr_test.pdf'),bbox_inches='tight');plt.close()

# # explore prior a bit
# count = prior.apply(func=lambda x:np.count_nonzero(x.values),axis=0,result_type='reduce')
# count.name = 'count'
# count = count.sort_values().to_frame().reset_index() # two column, index and count
# fig,ax = plt.subplots()
# ax.scatter(count.index.values,count['count'].values,s=2**2)
# texts = [ax.text(x=row.Index,y=row.count,s=row.index,fontsize=1,ha='center',va='center') for row in count.itertuples()]
# # adjust_text(texts)
# plt.savefig(os.path.join(root_path,'prior_scatter.pdf'),bbox_inches='tight');plt.close()

# sns.histplot(np.count_nonzero(prior.values,axis=0),bins=100)
# plt.savefig(os.path.join(root_path,'prior_test.pdf'),bbox_inches='tight');plt.close()

# start testing
# estimate = dc.run_aucell(mat=expr,net=net)
# estimate.to_csv(os.path.join(root_path,'decoupler','aucell.txt'),sep='\t')

# estimate, norm, pvals = dc.run_gsea(mat=expr,net=net,verbose=True)
# estimate.to_csv(os.path.join(root_path,'decoupler','gsea_estimate.txt'),sep='\t')
# norm.to_csv(os.path.join(root_path,'decoupler','gsea_norm.txt'),sep='\t')
# pvals.to_csv(os.path.join(root_path,'decoupler','gsea_pvals.txt'),sep='\t')

# estimate = dc.run_gsva(mat=expr,net=net,mx_diff=False)
# estimate.to_csv(os.path.join(root_path,'decoupler','gsva_es_method1.txt'),sep='\t')

# estimate = dc.run_mdt(mat=expr,net=net)
# estimate.to_csv(os.path.join(root_path,'mdt_estimate.txt'),sep='\t')

# estimate, pvals = dc.run_mlm(mat=expr,net=net)
# estimate.to_csv(os.path.join(root_path,'decoupler','mlm_estimate.txt'),sep='\t')

# estimate, pvals = dc.run_ora(mat=expr,net=net,n_background=35574)
# estimate.to_csv(os.path.join(root_path,'decoupler','ora_estimate.txt'),sep='\t')

# estimate = dc.run_udt(mat=expr,net=net)
# estimate.to_csv(os.path.join(root_path,'decoupler','udt_estimate.txt'),sep='\t')

# estimate, pvals = dc.run_ulm(mat=expr,net=net)
# estimate.to_csv(os.path.join(root_path,'decoupler','ulm_estimate.txt'),sep='\t')

# estimate, pvals = dc.run_viper(mat=expr,net=net)
# estimate.to_csv(os.path.join(root_path,'decoupler','viper_estimate.txt'),sep='\t')

# estimate, norm, corr, pvals = dc.run_wmean(mat=expr,net=net)
# estimate.to_csv(os.path.join(root_path,'decoupler','wmean_estimate.txt'),sep='\t')

# estimate, norm, corr, pvals = dc.run_wsum(mat=expr,net=net)
# estimate.to_csv(os.path.join(root_path,'decoupler','wsum_estimate.txt'),sep='\t')

# add metadata on MDT
sfa = pd.read_csv('mdt_estimate.txt',sep='\t',index_col=0).T
group = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/survival/groups.txt',sep='\t',index_col=0,header=None)
group = group[2].to_dict()  # 451 sample has burden label

burden = pd.read_csv('../TCGA_melanoma/result_new/burden_stage3.txt',sep='\t',index_col=0).iloc[-1,:][:-1]
burden = burden.to_dict()

sfa = sfa.loc[:,list(group.keys())] # 221 rows x 451 columns 

mi_array_sample = tuple(sfa.columns.tolist())
mi_array_group = tuple([group[item] for item in mi_array_sample])
mi_array_burden = tuple([burden[item] for item in mi_array_sample])
mi_array = [mi_array_sample,mi_array_group,mi_array_burden]
mi = pd.MultiIndex.from_arrays(arrays=mi_array,names=['sample','group','burden'])
sfa.columns = mi
# sfa.to_csv('mdt_estimate_morpheus.txt',sep='\t')

from scipy.stats import pearsonr, spearmanr, mannwhitneyu
import math
burden_label = mi.get_level_values(level='burden')
col1 = []
col2 = []
col3 = []
col4 = []
for i in range(sfa.shape[0]):
    s,p = pearsonr(burden_label,sfa.iloc[i,:].values)
    col1.append(s)
    col2.append(p)
    col3.append(-math.log10(p))
    h = sfa.iloc[i,].loc[(slice(None),'high',slice(None))]
    l = sfa.iloc[i,].loc[(slice(None),'low',slice(None))]
    s,p = mannwhitneyu(h,l)
    col4.append(p)
df = pd.DataFrame(index=sfa.index,data={'s':col1,'p':col2,'-log10':col3,'mann':col4}).sort_values(by='s').reset_index()
df.to_csv('result.txt',sep='\t')
df = pd.read_csv('result.txt',sep='\t',index_col=0)
df.sort_values(by='s',ascending=False,inplace=True)
fig,ax = plt.subplots(figsize=(10,4.8))
ax.scatter(x=df.index.tolist(),y=df['s'],s=1**2,c=['red' if item < 0.05 else 'black' for item in df['mann']])
ax.set_xticks(df.index.tolist())
ax.set_xticklabels(df['index'],fontsize=2,rotation=90)
ax.set_ylim([0,-0.35])
ax.set_xlabel('RBP')
ax.set_ylabel('negative correlation w/ burden')
plt.savefig('result.pdf',bbox_inches='tight')
plt.close()
sys.exit('stop')

'''
exp steady state                 Gene norm      
count steady state               Gene count
exp                              Junction norm
count                            Junction count

'''

# add exp
exp = pd.read_csv('/data/salomonis2/NCI-R01/TCGA-SKCM/TCGA_SKCM_hg38/ExpressionInput/exp.TCGA-SKCM-steady-state.txt',sep='\t',index_col=0)
ensg_dict = {
    'HNRNPA2B1':'ENSG00000122566',
    'HNRNPM':'ENSG00000099783',
    'SRSF1':'ENSG00000136450',
    'HNRNPK':'ENSG00000165119',
    'HNRNPA1':'ENSG00000135486',
    'U2AF1':'ENSG00000160201',
    'U2AF2':'ENSG00000063244',
    'SF3B1':'ENSG00000115524',
    'QKI':'ENSG00000112531',
    'SRSF5':'ENSG00000100650'
}
exp_subset = exp.loc[list(ensg_dict.values()),list(mi_array_sample)]
mi_exp = [tuple(exp_subset.iloc[i,:].values) for i in range(exp_subset.shape[0])]

# calculate association
for k,v in ensg_dict.items():
    s,p = pearsonr(exp_subset.loc[v,:].values,sfa.loc[k,:].values)
    print(k,v,s,p)

'''
HNRNPA2B1 ENSG00000122566 -0.3530657357102911 1.0974836902850253e-14
HNRNPM ENSG00000099783 -0.10864774850094357 0.02101258481101542
SRSF1 ENSG00000136450 -0.4756152698581512 7.780161279932694e-27
HNRNPK ENSG00000165119 0.003384684343442297 0.9428562079331373
HNRNPA1 ENSG00000135486 -0.14224685563981398 0.0024626771659759525
U2AF1 ENSG00000160201 -0.09302425118554272 0.04834375633303713
U2AF2 ENSG00000063244 0.02900441583527343 0.5389644364554604
SF3B1 ENSG00000115524 -0.507938216846988 5.92314283252423e-31
QKI ENSG00000112531 -0.17983195144079786 0.00012313570313980186
SRSF5 ENSG00000100650 -0.27794654129833923 1.9129526235481513e-09
'''


# add mutation
mut = pd.read_csv('../TCGA_melanoma/TCGA-SKCM.mutect2_snv.tsv',sep='\t',index_col=0)
mut_sub = mut.loc[mut['gene']=='SF3B1',:]
all_samples = set(mut_sub.index)
mi_mut = tuple([True if '-'.join(item.split('-')[:4]) in all_samples else False for item in mi_array_sample])

mi = mi_exp + [mi_array_sample,mi_array_group,mi_array_burden] + [mi_mut]
mi = pd.MultiIndex.from_arrays(arrays=mi,names=list(ensg_dict.keys()) + ['sample','group','burden'] + ['mut'])
sfa.columns = mi
sfa.to_csv('mdt_estimate_morpheus_add_gene_mutation.txt',sep='\t')



























































