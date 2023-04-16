#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os,subprocess,sys
from scipy.stats import norm
from tqdm import tqdm
import pickle
import seaborn as sns
import matplotlib.pyplot as plt




def get_sfa_and_full_prior(expr,prior):
    n_event, n_sample = expr.shape
    n_event_subset, n_sf = prior.shape
    extra_events = list(set(expr.index).difference(set(prior.index)))  # n_event-n_event_subset
    enlarge_prior = pd.DataFrame(data=np.zeros((len(extra_events),n_sf)),columns=prior.columns,index=extra_events) # n_event-n_event_subset, n_sf
    full_prior = pd.concat([prior,enlarge_prior],axis=0).loc[expr.index,:]  # n_event, n_sf
    sfa_data = np.linalg.pinv(full_prior.values)@expr.values
    # sfa_data = (sfa_data - sfa_data.mean(axis=1)[:,np.newaxis]) / sfa_data.std(axis=1)[:,np.newaxis]   # scale
    sfa = pd.DataFrame(data=sfa_data,columns=expr.columns,index=full_prior.columns) # n_sf * n_sample
    # sfa = sfa.apply(func=lambda x:(x - x.min())/(x.max() - x.min()),axis=1,result_type='expand')
    # sns.clustermap(sfa,method='average',metric='correlation',cmap='viridis')
    # plt.savefig('../altanalyze_output/inference/bbsr_weight1/sfa_clustermap_{}.pdf'.format(cutoff),bbox_inches='tight')
    # plt.close()
    assert np.all(expr.index == full_prior.index)
    return sfa, full_prior


def SSE_pred(X_test,X_train,B_train,A_test,A_train,cutoff=None,normalize_B=False):
    '''
    X_test: n_event * n_sample_test
    X_train: n_event * n_sample_train
    B_train: n_event * n_sf      # n_nonzero_samples
    A_test: n_sf * n_sample_test
    A_train: n_sf * n_sample_train
    '''
    if cutoff is not None:
        B_train.where(B_train > cutoff,other=0,inplace=True)
    if normalize_B:
        max_ = np.amax(B_train.values)
        min_ = np.amin(B_train.values)
        B_train = pd.DataFrame(data=(B_train.values - min_)/(max_-min_),columns=B_train.columns,index=B_train.index)
    E = X_test.values - B_train.values @ ((A_test.values - A_train.values.mean(axis=1)[:,np.newaxis])/A_train.values.std(axis=1)[:,np.newaxis]) - X_train.values.mean(axis=1)[:,np.newaxis]
    # E = X_test.values - B_train.values @ A_test.values
    SE = E ** 2
    SSE = SE.sum()
    return SSE

def SSE_null(X_test,X_train):
    E = X_test.values - X_train.values.mean(axis=1)[:,np.newaxis]
    # E = X_test.values - X_test.values.mean(axis=1)[:,np.newaxis]
    SE = E ** 2
    SSE = SE.sum()
    return SSE

def R2(pred,null):
    r2 = 1 - pred/null
    return r2

def run(cutoff):
    print(cutoff)
    A_train,_ = get_sfa_and_full_prior(train,prior)
    A_test,_ = get_sfa_and_full_prior(test,prior)

    sse_pred = SSE_pred(test,train,network,A_test,A_train,cutoff=cutoff)
    print('SSE_pred:',sse_pred)
    sse_null = SSE_null(test,train)
    print('SSE_null:',sse_null)
    r2 = R2(sse_pred,sse_null)
    print('R2:',r2)


train = pd.read_csv('../altanalyze_output/inference/expr_train.tsv',sep='\t',index_col=0).T
test = pd.read_csv('../altanalyze_output/inference/expr_test.tsv',sep='\t',index_col=0).T
network = pd.read_csv('../altanalyze_output/inference/bbsr_weight1/network_reformat.txt',sep='\t',index_col=0)
prior = pd.read_csv('../altanalyze_output/inference/prior.tsv',sep='\t',index_col=0)

for cutoff in [0,0.2,0.4,0.6,0.8,0.9]:
    run(cutoff)
























