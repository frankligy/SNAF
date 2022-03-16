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


def gtex_configuration(gtex_db,t_min_arg,n_max_arg,add_control=None):
    global adata
    global adata_gtex
    global t_min
    global n_max
    adata_gtex = anndata.read_h5ad(gtex_db)
    if add_control is not None:
        print('adding additional control samples {} to the database'.format(add_control.shape))
        tissue_dict = adata_gtex.var['tissue'].to_dict()
        tissue_dict_right = {k:'additional_control' for k in add_control.columns}
        tissue_dict.update(tissue_dict_right)
        df_left = adata_gtex.to_df()
        df_right = add_control
        df_combine = df_left.join(other=df_right,how='outer').fillna(0)
        adata = anndata.AnnData(X=df_combine.values,obs=pd.DataFrame(index=df_combine.index),var=pd.DataFrame(index=df_combine.columns))
        print('now the shape of control db is {}'.format(adata.shape))
        adata.var['tissue'] = adata.var_names.map(tissue_dict).values
        adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
        total_count = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
        adata.var['total_count'] = total_count
    else:
        adata = adata_gtex
    t_min = t_min_arg
    n_max = n_max_arg




def multiple_crude_sifting(junction_count_matrix,add_control=None):   # for JunctionCountMatrixQuery class, only consider gtex
    df = pd.DataFrame(index=junction_count_matrix.index,data = {'max':junction_count_matrix.max(axis=1).values})
    # consider gtex
    junction_to_mean = adata_gtex.obs.loc[adata_gtex.obs_names.isin(junction_count_matrix.index),'mean'].to_dict()
    df['mean'] = df.index.map(junction_to_mean).fillna(value=0)
    df['diff'] = df['max'] - df['mean']
    df['cond'] = (df['mean'] < n_max) & (df['diff'] > t_min)
    valid = df.loc[df['cond']].index.tolist()
    # consider add_control
    if add_control is not None:
        junction_to_mean = add_control.mean(axis=1).to_dict()
        df['mean_add'] = df.index.map(junction_to_mean).fillna(value=0)
        df['diff_add'] = df['max'] - df['mean_add']
        df['cond_add'] = (df['mean_add'] < n_max) & (df['diff_add'] > t_min)
        valid_add = df.loc[df['cond_add']].index.tolist()
        valid = list(set(valid).intersection(set(valid_add)))
    invalid = list(set(junction_count_matrix.index).difference(set(valid)))
    # now, consider each entry
    gtex_df = pd.concat([df['mean']]*junction_count_matrix.shape[1],axis=1)
    gtex_df.columns = junction_count_matrix.columns
    diff_df_gtex = junction_count_matrix - gtex_df
    cond_df = (gtex_df < n_max) & (diff_df_gtex > t_min)
    if add_control is not None:
        add_df = pd.concat([df['mean_add']]*junction_count_matrix.shape[1],axis=1)    
        add_df.columns = junction_count_matrix.columns
        diff_df_add = junction_count_matrix - add_df
        cond_df = cond_df & (add_df < n_max) & (diff_df_add > t_min)
    return valid,invalid,cond_df


def crude_tumor_specificity(uid,count):    # for NeoJunction class, since we normally start from Jcmq with check_gtex=False, rarely being called.
    detail = ''
    if uid not in set(adata_gtex.obs_names):
        mean_value = 0
    else:
        mean_value = adata_gtex.obs.loc[uid,'mean']
    diff = count - mean_value
    if mean_value < n_max and diff >= t_min:
        identity = True
    else:
        identity = False
    return identity,mean_value


def mle_func(parameters,y):
    sigma = parameters
    ll = np.sum(stats.halfnorm.logpdf(y,0,sigma))
    neg_ll = -1 * ll
    return neg_ll

def tumor_specificity(uid,method,return_df=False):
    try:
        info = adata[[uid],:]
    except:
        print('{} not detected in gtex, impute as zero'.format(query))
        info_tmp = adata[['ENSG00000090339:E4.3-E4.5'],:]
        info = ad.AnnData(X=csr_matrix(np.full((1,info_tmp.shape[1]),0)),obs=info_tmp.obs,var=info_tmp.var)  # weired , anndata 0.7.6 can not modify the X in place? anndata 0.7.2 can do that in scTriangulate
    df = pd.DataFrame(data={'value':info.X.toarray().squeeze(),'tissue':info.var['tissue'].values},index=info.var_names)
    if method == 'mean':
        try:
            sigma = adata.obs.loc[uid,'mean']
        except KeyError:
            sigma = 0
        if return_df:
            return sigma,df
        else:
            return sigma
    elif method == 'mle':
        try:
            y = adata[[uid],:].X.toarray().squeeze() / adata.var['total_count'].values
        except KeyError:
            sigma = 0
        else:
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
        if return_df:
            return sigma,df
        else:
            return sigma
    elif method == 'bayesian':
        try:
            y = adata[[uid],:].X.toarray().squeeze() / adata.var['total_count'].values
        except KeyError:
            sigma = 0
        else:
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
                trace = pm.sample(100,step=step,return_inferencedata=False,cores=1)
            df = az.summary(trace,round_to=2)
            sigma = df.iloc[0]['mean']
        if return_df:
            return sigma,df
        else:
            return sigma
            
















