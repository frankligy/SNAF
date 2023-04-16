#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os,subprocess,sys
from scipy.stats import norm
from tqdm import tqdm
import pickle
import plotly.graph_objects as go
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt

# ## Step1: no threshold, how large is the network
# network = pd.read_csv('../altanalyze_output/HepG2_inference/bbsr_weight1/network_aug.txt',sep='\t',index_col=0)
# network = network.stack().reset_index()
# network.columns = ['event','sf','weight']
# network = network.loc[network['weight']>0,:]
# network.sort_values(by='weight',ascending=False,inplace=True)

# # look from SF perspective
# size_dict = {}
# for sf, sub_df in network.groupby(by='sf'):
#     size_dict[sf] = sub_df.shape[0]
# size_dict = {k:v for k,v in sorted(size_dict.items(),key=lambda x:x[1])}

# value_dict = {}
# for sf, sub_df in network.groupby(by='sf'):
#     value_dict[sf] = np.median(sub_df['weight'].values)
# value_dict = {k:v for k,v in sorted(value_dict.items(),key=lambda x:x[1])}

# order = list(value_dict.keys())
# fig,ax = plt.subplots(figsize=(6.4,30))
# sns.violinplot(x='weight',y='sf',data=network,order=order,ax=ax,size=0.5)
# ylabel = ax.yaxis.get_ticklabels()
# ylabel = [text.get_text() + '_' + str(size_dict[text.get_text()]) for text in ylabel]
# ax.set_yticklabels(ylabel)
# plt.savefig('../altanalyze_output/HepG2_inference/bbsr_weight1/threshold/sf_value.pdf',bbox_inches='tight')
# plt.close()


# # look from event perspective
# size_dict = {}
# for event, sub_df in network.groupby(by='event'):
#     size_dict[event] = sub_df.shape[0]
# size_dict = {k:v for k,v in sorted(size_dict.items(),key=lambda x:x[1])}
# se = pd.Series(data=size_dict,name='size')
# fig,ax = plt.subplots()
# ax.hist(se.values,bins=50,edgecolor='k')
# plt.savefig('../altanalyze_output/HepG2_inference/bbsr_weight1/threshold/event_size.pdf',bbox_inches='tight')
# plt.close()



### Step2: implement MCC and F1, precision and recell, etc
network = pd.read_csv('../altanalyze_output/HepG2_inference/bbsr_weight1/network_aug.txt',sep='\t',index_col=0)
network = network.stack().to_frame(name='network')
network.index = [event + ',' + sf for event,sf in network.index.tolist()]

prior = pd.read_csv('../altanalyze_output/HepG2_inference/prior.tsv',sep='\t',index_col=0)
prior = prior.stack().to_frame(name='prior')
prior.index = [event + ',' + sf for event,sf in prior.index.tolist()]

def thresholding(network,prior,interval,method,outdir):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    from sklearn.metrics import matthews_corrcoef, f1_score, precision_score, recall_score
    merged = network.join(prior,how='left')
    merged.sort_values(by='network',inplace=True,ascending=False)
    merged = merged.loc[merged['network']>0,:]
    merged['prior'].fillna(value=0,inplace=True)
    merged['prior'].where(cond=merged['prior']==0,other=1,inplace=True)
    conf_list, metric_list = [], []
    for i in tqdm(range(0,merged.shape[0],interval),total=merged.shape[0]//interval):
        conf = merged['network'].iloc[i]
        merged_copy = merged.copy()
        merged_copy['network'].where(cond=merged_copy['network']>=conf,other=0,inplace=True)
        merged_copy['network'].where(cond=merged_copy['network']<conf,other=1,inplace=True)
        if method == 'mcc':
            metric = matthews_corrcoef(merged_copy['prior'].values,merged_copy['network'].values)
        elif method == 'f1':
            metric = f1_score(merged_copy['prior'].values,merged_copy['network'].values)
        elif method == 'precision':
            metric = precision_score(merged_copy['prior'].values,merged_copy['network'].values)
        elif method == 'recall':
            metric = recall_score(merged_copy['prior'].values,merged_copy['network'].values)
        conf_list.append(conf)
        metric_list.append(metric)
    fig,ax = plt.subplots()
    ax.plot(conf_list,metric_list)
    ax.set_xlim([1,0])
    plt.savefig(os.path.join(outdir,'threshold_{}.pdf'.format(method)),bbox_inches='tight')
    plt.close()

thresholding(network,prior,10000,'recall',outdir='../altanalyze_output/HepG2_inference/bbsr_weight1/threshold')