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


# # test on HepG2
# # expr: sample * event
# # prior: event * sf
expr_path = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/HepG2_inference/expr.tsv'
prior_path = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/HepG2_inference/prior.tsv'
root_path = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/HepG2_inference'
expr = pd.read_csv(expr_path,sep='\t',index_col=0)  # 446 rows x 35573 columns
prior = pd.read_csv(prior_path,sep='\t',index_col=0)  # 73157 rows x 221 columns
common_events = list(set(expr.columns).intersection(set(prior.index)))
expr = expr.loc[:,common_events]   # 446 rows x 28733 columns
prior = prior.loc[common_events,:]  # 28733 rows x 221 columns
net = prior.stack().reset_index()
net.columns = ['target','source','weight']
net = net.loc[net['weight']!=0,:]

# explore expr a bit
# sns.histplot(np.count_nonzero(expr.values,axis=1),bins=100)
# plt.savefig(os.path.join(root_path,'expr_test.pdf'),bbox_inches='tight');plt.close()

# explore prior a bit
count = prior.apply(func=lambda x:np.count_nonzero(x.values),axis=0,result_type='reduce')
count.name = 'count'
count = count.sort_values().to_frame().reset_index() # two column, index and count
fig,ax = plt.subplots()
ax.scatter(count.index.values,count['count'].values,s=2**2)
texts = [ax.text(x=row.Index,y=row.count,s=row.index,fontsize=1,ha='center',va='center') for row in count.itertuples()]
# adjust_text(texts)
plt.savefig(os.path.join(root_path,'prior_scatter.pdf'),bbox_inches='tight');plt.close();sys.exit('stop')

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
# estimate.to_csv(os.path.join(root_path,'decoupler','mdt_estimate.txt'),sep='\t')

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

# evaluate
def add_meta_data(name,sfa=None,actual_name=None):
    if name is None:  # for bbsr and stars
        sfa = sfa
        name = actual_name
    else:
        sfa = pd.read_csv(os.path.join(root_path,'decoupler','{}.txt'.format(name)),sep='\t',index_col=0).T
    meta_k = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
    meta_c = pd.read_csv('../processed_metadata_control.tsv',sep='\t',index_col=0)
    meta = pd.concat([meta_k,meta_c])
    meta.columns = [item.replace(' ','_') for item in meta.columns]
    meta['Experiment_target'].fillna(value='control',inplace=True)
    meta['Experiment_target'] = [item.split('-')[0] if item != 'control' else item for item in meta['Experiment_target']]
    pair1_to_target = meta['Experiment_target'].to_dict()
    mi_array_target = tuple([pair1_to_target[item.split('_')[0]] for item in sfa.columns])
    mi_array_name = tuple(sfa.columns.tolist())
    rep_counter = {}
    rep = []
    for t,n in zip(mi_array_target,mi_array_name):
        if t == 'control':
            rep.append('c')
        else:
            if rep_counter.get(t,None) is None:
                rep_counter[t] = 1
            else:
                rep_counter[t] += 1
            count = rep_counter[t]
            rep.append('rep{}'.format(count))
    mi_array_rep = rep
    mi = pd.MultiIndex.from_arrays(arrays=[mi_array_rep,mi_array_target,mi_array_name],names=('rep','target','name'))
    sfa.columns = mi
    sfa.to_csv(os.path.join(root_path,'decoupler','{}_visual_morpheus.txt'.format(name)),sep='\t')
    return sfa

def calculate_metrics(sfa,name):
    sfa = sfa.copy()
    common_sf = list(set(sfa.columns.get_level_values('target')).intersection(sfa.index))
    sfa_quantify = sfa.loc[common_sf,(slice(None),common_sf+['control'],slice(None))]
    from sklearn.metrics import coverage_error, label_ranking_average_precision_score, label_ranking_loss, ndcg_score
    score_dict = {'ce':0,'lraps':0,'lrl':0,'ndcgs':0}
    score_dict_sf = {}
    for sf in tqdm(sfa_quantify.index,total=sfa_quantify.shape[0]):
        tmp = sfa_quantify.loc[[sf],(slice(None),[sf]+['control'],slice(None))]
        tmp.columns = tmp.columns.get_level_values('target')
        y_score = np.nan_to_num(np.negative(tmp.values),nan=0,posinf=1000,neginf=-1000)   # there are some score is labelled as inf
        y_true = tmp.columns.map({sf:1,'control':0}).values.reshape(1,-1)
        all_metrics = []
        # compute coverage_error
        metric = coverage_error(y_true,y_score)
        score_dict['ce'] += metric
        all_metrics.append(metric)
        # compute label ranking average precision score
        metric = label_ranking_average_precision_score(y_true,y_score)
        score_dict['lraps'] += metric
        all_metrics.append(metric)
        # compute label ranking loss
        metric = label_ranking_loss(y_true,y_score)
        score_dict['lrl'] += metric
        all_metrics.append(metric)
        # compute normalized discounted cumulative gain
        metric = ndcg_score(y_true,y_score)
        score_dict['ndcgs'] += metric
        all_metrics.append(metric)
        # summary
        score_dict_sf[sf] = all_metrics
    score_dict['ce'] = score_dict['ce']/sfa_quantify.shape[0]
    score_dict['lraps'] = score_dict['lraps']/sfa_quantify.shape[0]
    score_dict['lrl'] = score_dict['lrl']/sfa_quantify.shape[0]
    score_dict['ndcgs'] = score_dict['ndcgs']/sfa_quantify.shape[0]
    df = pd.DataFrame(data=score_dict_sf,index=['coverage_error','label_ranking_average_precision_score','label_ranking_loss','normalized_discounted_cumulative_gain']).T
    df.to_csv(os.path.join(root_path,'decoupler','{}_evaluation.txt'.format(name)),sep='\t')
    return df

methods = ['aucell','viper_estimate','mdt_estimate','wsum_estimate','wmean_estimate','ulm_estimate','ora_estimate','mlm_estimate','gsva_es_method1','gsva_es_method2','gsea_estimate','gsea_norm']
# for m in methods:
#     sfa = add_meta_data(m)
#     df = calculate_metrics(sfa,m)

# sfa = pd.read_csv('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/HepG2_inference/stars_b0.5/sfa_network.txt',index_col=0,sep='\t',skiprows=1,header=0)
# sfa.index.name = None
# m = 'star_0.5'
# sfa = add_meta_data(None,sfa,m)
# df = calculate_metrics(sfa,m)

# sfa = pd.read_csv('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/HepG2_inference/bbsr_weight1/sfa_network.txt',index_col=0,sep='\t',skiprows=1,header=0)
# sfa.index.name = None
# m = 'bbsr_1'
# sfa = add_meta_data(None,sfa,m)
# df = calculate_metrics(sfa,m)

# combine evaluation
methods_plus = methods + ['star_0.5','bbsr_1']
# data = []
# for m in methods_plus:
#     df = pd.read_csv(os.path.join(root_path,'decoupler','{}_evaluation.txt'.format(m)),sep='\t',index_col=0)
#     for row in df.itertuples():
#         data.append((row.Index,'coverage_error',row.coverage_error,m))
#         data.append((row.Index,'label_ranking_average_precision_score',row.label_ranking_average_precision_score,m))
#         data.append((row.Index,'label_ranking_loss',row.label_ranking_loss,m))
#         data.append((row.Index,'normalized_discounted_cumulative_gain',row.normalized_discounted_cumulative_gain,m))
# combine = pd.DataFrame.from_records(data=data,columns=['sf','metric','value','method'])
# combine.to_csv(os.path.join(root_path,'decoupler','combined_evaluation.txt'),sep='\t')

# visualize one, boxplot
# combine_ori = pd.read_csv(os.path.join(root_path,'decoupler','combined_evaluation.txt'),sep='\t',index_col=0)
# sfs = ['SF3B1','U2AF2','U2AF1','FUS','HNRNPK']
# combine_ori = combine_ori.loc[combine_ori['sf'].isin(sfs),:]
# for metric in ['coverage_error','label_ranking_average_precision_score','label_ranking_loss','normalized_discounted_cumulative_gain']:
#     combine = combine_ori.loc[combine_ori['metric']==metric,:].copy()
#     custom_order = combine.groupby(by='method').apply(lambda x:x['value'].mean()).sort_values(ascending=False).index.tolist()
#     combine['method'] = combine['method'].astype(dtype=pd.CategoricalDtype(categories=custom_order,ordered=True))
#     combine.sort_values(by='method',inplace=True)
#     fig,ax = plt.subplots()
#     sns.boxplot(data=combine,x='method',y='value',ax=ax)
#     ax.tick_params(labelsize=2)
#     ax.set_title(metric)
#     plt.savefig('test_{}.pdf'.format(metric),bbox_inches='tight');plt.close()

# visualize two, heatmap
metric = 'label_ranking_average_precision_score'
lis = []
for m in methods_plus:
    df = pd.read_csv(os.path.join(root_path,'decoupler','{}_evaluation.txt'.format(m)),sep='\t',index_col=0)
    lis.append(df[metric])
combine = pd.concat(lis,axis=1)
combine.columns = methods_plus
sort_sf = combine.apply(lambda x:x.mean(),axis=1).sort_values(ascending=False).index
combine = combine.loc[sort_sf,:]
combine.to_csv(os.path.join(root_path,'decoupler','heatmap_{}_all_sfs.txt'.format(metric)),sep='\t')
























































