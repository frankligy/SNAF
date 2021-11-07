#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt
import anndata
from scipy.optimize import minimize
from scipy import stats
import pymc3 as pm   # conda install -c conda-forge pymc3 mkl-service
import theano
import arviz as az

'''
this script is to query the tumor specificity of the junction
'''


def gtex_configuration(gtex_db,t_min_arg,n_max_arg):
    global adata
    global t_min
    global n_max
    adata = anndata.read_h5ad(gtex_db)
    t_min = t_min_arg
    n_max = n_max_arg


def multiple_crude_sifting(junction_count_matrix):
    df = pd.DataFrame(index=junction_count_matrix.index,
                      data = {'max':junction_count_matrix.max(axis=1).values})
    junction_to_mean = adata.obs.loc[adata.obs_names.isin(junction_count_matrix.index),'mean'].to_dict()
    df['mean'] = df.index.map(junction_to_mean).fillna(value=0)
    df['diff'] = df['max'] - df['mean']
    df['cond'] = (df['mean'] < n_max) & (df['diff'] > t_min)
    valid = df.loc[df['cond']].index.tolist()
    invalid = df.loc[~df['cond']].index.tolist()

    gtex_df = pd.concat([df['mean']]*junction_count_matrix.shape[1],axis=1)
    gtex_df.columns = junction_count_matrix.columns
    diff_df = junction_count_matrix - gtex_df
    cond_df = (gtex_df < n_max) & (diff_df > t_min)
    return valid,invalid,cond_df





def crude_tumor_specificity(uid,count):
    if uid not in set(adata.obs_names):
        mean_value = 0
    else:
        mean_value = adata.obs.loc[uid,'mean']
    diff = count - mean_value
    if mean_value < n_max and diff >= t_min:
        identity = True
    else:
        identity = False
    return identity


def mle_func(parameters,y):
    sigma = parameters
    ll = np.sum(stats.halfnorm.logpdf(y,0,sigma))
    neg_ll = -1 * ll
    return neg_ll

def accurate_tumor_specificity(uid,method):
    adata = anndata.read_h5ad('../data/GTEx_junction_counts.h5ad')
    if method == 'mle':
        y = adata[[uid],:].X.toarray().squeeze() / adata.var['total_count'].values
        mle_model = minimize(mle_func,np.array([0.2]),args=(y,),bounds=((0,1),),method='L-BFGS-B')
        '''
        fun: nan
        hess_inv: <1x1 LbfgsInvHessProduct with dtype=float64>
        jac: array([6368.36862213])
        message: b'ABNORMAL_TERMINATION_IN_LNSRCH'
        nfev: 42
        nit: 0
        status: 2
        success: False
        x: array([0.2])
        '''
        if mle_model.success:
            sigma = mle_model.x[0]
        else:   # usually means too many zero, so true expression is near zero
            sigma = 0
    elif method == 'bayesian':
        y = adata[[uid],:].X.toarray().squeeze() / adata.var['total_count'].values
        x = []
        for tissue in adata.var['tissue'].unique():
            sub = adata[uid,adata.var['tissue']==tissue]
            c = np.count_nonzero(sub.X.toarray())
            x.append(c)
        x = np.array(x)
        with pm.Model() as m:
            sigma = pm.Uniform('sigma',lower=0,upper=1)
            nc = pm.HalfNormal('nc',sigma=sigma,observed=y)
            rate = pm.math.sum(nc)/len(y)
            c = pm.Poisson('c',mu=rate,observed=x)
        with m:
            step = pm.NUTS()
            trace = pm.sample(200,step=step,return_inferencedata=False,cores=1)
        df = az.summary(trace,round_to=2)
        sigma = df.iloc[0]['mean']
    return 1-sigma
            
















