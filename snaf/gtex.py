#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt
import anndata as ad
from scipy.optimize import minimize, minimize_scalar
from scipy import stats
from scipy.sparse import csr_matrix, find
from tqdm import tqdm
import re

try:
    import pymc3 as pm   # conda install -c conda-forge pymc3 mkl-service
    import theano
    import arviz as az
except ImportError:
    print('''
        Optional package pymc3 is not installed, it is for calculating tumor specificity using hirerarchical bayesian model
        For Linux: https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(Linux)
        For MacOS: https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(MacOS)
        For PC:    https://github.com/pymc-devs/pymc/wiki/Installation-Guide-(Windows)
    ''')

'''
this script is to query the tumor specificity of the junction
'''


def gtex_configuration(df,gtex_db,t_min_arg,n_max_arg,normal_cutoff_arg,tumor_cutoff_arg,normal_prevalance_cutoff_arg,tumor_prevalance_cutoff_arg,add_control=None):
    global adata_gtex
    global adata
    global t_min
    global n_max
    global normal_cutoff
    global tumor_cutoff
    global normal_prevalance_cutoff
    global tumor_prevalance_cutoff
    tested_junctions = set(df.index)
    adata = ad.read_h5ad(gtex_db)
    adata = adata[np.logical_not(adata.obs_names.duplicated()),:] 
    adata = adata[list(set(adata.obs_names).intersection(tested_junctions)),:]  
    print('Current loaded gtex cohort with shape {}'.format(adata.shape))
    tissue_dict = adata.var['tissue'].to_dict()
    adata_gtex = adata   # already has mean and tissue variables
    if add_control is not None:
        for id_, control in add_control.items():
            if isinstance(control,pd.DataFrame):
                assert len(set(control.columns).intersection(tissue_dict.keys())) == 0  # sample id can not be ambiguous
                control = control.loc[np.logical_not(control.index.duplicated()),:]
                control = control.loc[list(set(control.index).intersection(tested_junctions)),:]
                print('Adding cohort {} with shape {} to the database'.format(id_,control.shape))
                tissue_dict_right = {k:id_ for k in control.columns}
                tissue_dict.update(tissue_dict_right)
                df_left = adata.to_df()
                df_right = control
                df_combine = df_left.join(other=df_right,how='outer').fillna(0)
                adata = ad.AnnData(X=csr_matrix(df_combine.values),obs=pd.DataFrame(index=df_combine.index),var=pd.DataFrame(index=df_combine.columns))

            elif isinstance(control,ad.AnnData):
                assert len(set(control.var_names).intersection(tissue_dict.keys())) == 0
                control = control[np.logical_not(control.obs_names.duplicated()),:]
                control = control[list(set(control.obs_names).intersection(tested_junctions)),:]
                print('Adding cohort {} with shape {} to the database'.format(id_,control.shape))
                if 'tissue' in control.var.columns:   # if tissue is in var columns, it will be used 
                    tissue_dict_right = control.var['tissue'].to_dict()
                else:
                    tissue_dict_right = {k:id_ for k in control.var_names}
                tissue_dict.update(tissue_dict_right)
                df_left = adata.to_df()
                df_right = control.to_df()
                df_combine = df_left.join(other=df_right,how='outer').fillna(0)
                adata = ad.AnnData(X=csr_matrix(df_combine.values),obs=pd.DataFrame(index=df_combine.index),var=pd.DataFrame(index=df_combine.columns))
            
            else:
                raise Exception('control must be either in dataframe or anndata format')

            print('now the shape of control db is {}'.format(adata.shape))
            adata.var['tissue'] = adata.var_names.map(tissue_dict).values
            adata.obs['mean'] = np.array(adata.X.mean(axis=1)).squeeze()
            total_count = np.array(adata.X.sum(axis=0)).squeeze() / 1e6
            adata.var['total_count'] = total_count
            

    t_min = t_min_arg
    n_max = n_max_arg
    normal_cutoff = normal_cutoff_arg
    tumor_cutoff = tumor_cutoff_arg
    normal_prevalance_cutoff = normal_prevalance_cutoff_arg
    tumor_prevalance_cutoff = tumor_prevalance_cutoff_arg

    return adata


def multiple_crude_sifting(junction_count_matrix,add_control,dict_exonlist,outdir,filter_mode):
    if filter_mode == 'prevalance':
        valid,invalid,cond_df = multiple_crude_sifting_prevalance(junction_count_matrix,add_control,dict_exonlist,outdir)
    elif filter_mode == 'maxmin':
        valid,invalid,cond_df = multiple_crude_sifting_maxmin(junction_count_matrix,add_control,dict_exonlist,outdir)
    return valid,invalid, cond_df

def multiple_crude_sifting_prevalance(junction_count_matrix,add_control=None,dict_exonlist=None,outdir='.'):

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    df_to_write = []
    df = pd.DataFrame(index=junction_count_matrix.index)
    prevalance_tumor = np.count_nonzero((junction_count_matrix > tumor_cutoff).values,axis=1) / junction_count_matrix.shape[1]
    df['prevalance_tumor'] = prevalance_tumor
    # consider gtex
    prevalance_normal = np.count_nonzero((adata_gtex.X > normal_cutoff).toarray(),axis=1) / adata_gtex.shape[1]
    prevalance_normal_dict = {j:v for j,v in zip(adata_gtex.obs_names,prevalance_normal)}
    df['prevalance_normal'] = df.index.map(prevalance_normal_dict).fillna(value=0)
    df['cond'] = (df['prevalance_tumor'] > tumor_prevalance_cutoff) & (df['prevalance_normal'] < normal_prevalance_cutoff)
    valid = df.loc[df['cond']].index.tolist()
    tmp = df.copy()
    df_to_write.append(tmp)
    print('reduce valid NeoJunction from {} to {} because they are present in GTEx'.format(df.shape[0],len(valid)))
    if dict_exonlist is not None:   # a valid junction can not be present in any ensembl documented transcript
        updated_valid = []
        for uid in tqdm(valid):
            ensg = uid.split(':')[0]
            exons = ':'.join(uid.split(':')[1:])
            if '_' in exons or 'U' in exons or 'ENSG' in exons or 'I' in exons:
                updated_valid.append(uid)
            else:
                exonlist = dict_exonlist[ensg]
                exonstring = '|'.join(exonlist)
                e1,e2 = exons.split('-')
                pattern1 = re.compile(r'^{}\|{}\|'.format(e1,e2))  # ^E1.1|E2.3|
                pattern2 = re.compile(r'\|{}\|{}$'.format(e1,e2))  # |E1.1|E2.3$
                pattern3 = re.compile(r'\|{}\|{}\|'.format(e1,e2)) # |E1.1|E2.3|
                if re.search(pattern3,exonstring) or re.search(pattern2,exonstring) or re.search(pattern1,exonstring):   # as long as match one pattern, should be eliminated
                    continue
                else:
                    updated_valid.append(uid)
        print('reduce valid Neojunction from {} to {} because they are present in Ensembl db'.format(len(valid),len(updated_valid)))
        valid = updated_valid
    # consider add_control
    if add_control is not None:
        for i,(id_,control) in enumerate(add_control.items()):
            n_previous_valid = len(valid)
            if isinstance(control,pd.DataFrame):
                prevalance_normal = np.count_nonzero((control > normal_cutoff).values,axis=1) / control.shape[1]
                prevalance_normal_dict = {j:v for j,v in zip(control.index, prevalance_normal)}
            elif isinstance(control,ad.AnnData):
                prevalance_normal = np.count_nonzero((control.X > normal_cutoff).toarray(),axis=1) / control.shape[1]
                prevalance_normal_dict = {j:v for j,v in zip(control.obs_names,prevalance_normal)}
            else:
                raise Exception('control must be either in dataframe or anndata format')
            df['prevalance_normal_add'] = df.index.map(prevalance_normal_dict).fillna(value=0)
            df['cond_add'] = (df['prevalance_tumor'] > tumor_prevalance_cutoff) & (df['prevalance_normal_add'] < normal_prevalance_cutoff)
            valid_add = df.loc[df['cond_add']].index.tolist()
            valid = list(set(valid).intersection(set(valid_add)))
            tmp = df.copy(); tmp.drop(columns=['prevalance_tumor','prevalance_normal','cond'],inplace=True); tmp.rename(columns=lambda x:x+'_{}'.format(id_),inplace=True)
            df_to_write.append(tmp)
            print('reduce valid Neojunction from {} to {} because they are present in added control {}'.format(n_previous_valid,len(valid),id_))
    invalid = list(set(junction_count_matrix.index).difference(set(valid)))
    # now, consider each entry
    t_min = tumor_cutoff
    valid_set = set(valid)
    cond_dict = {j:(True if j in valid_set else False) for j in junction_count_matrix.index}
    tmp = pd.DataFrame(index=junction_count_matrix.index,data={'placeholder':junction_count_matrix.index.map(cond_dict).values})
    first_half_cond_df = pd.concat([tmp]*junction_count_matrix.shape[1],axis=1)
    first_half_cond_df.columns = junction_count_matrix.columns
    cond_df = (first_half_cond_df) & (junction_count_matrix > t_min)
    # write the df
    df_to_write = pd.concat(df_to_write,axis=1)
    df_to_write.to_csv(os.path.join(outdir,'NeoJunction_statistics_prevalance.txt'),sep='\t')
    return valid,invalid,cond_df



def multiple_crude_sifting_maxmin(junction_count_matrix,add_control=None,dict_exonlist=None,outdir='.'):   # for JunctionCountMatrixQuery class, only consider gtex
    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    df = pd.DataFrame(index=junction_count_matrix.index,data = {'max':junction_count_matrix.max(axis=1).values})
    df_to_write = []
    # consider gtex
    junction_to_mean = adata_gtex.obs.loc[adata_gtex.obs_names.isin(junction_count_matrix.index),'mean'].to_dict()
    df['mean'] = df.index.map(junction_to_mean).fillna(value=0)
    df['diff'] = df['max'] - df['mean']
    df['cond'] = (df['mean'] < n_max) & (df['diff'] > t_min)
    valid = df.loc[df['cond']].index.tolist()
    tmp = df.copy()
    df_to_write.append(tmp)
    print('reduce valid NeoJunction from {} to {} because they are present in GTEx'.format(df.shape[0],len(valid)))
    if dict_exonlist is not None:   # a valid junction can not be present in any ensembl documented transcript
        updated_valid = []
        for uid in tqdm(valid):
            ensg = uid.split(':')[0]
            exons = ':'.join(uid.split(':')[1:])
            if '_' in exons or 'U' in exons or 'ENSG' in exons or 'I' in exons:
                updated_valid.append(uid)
            else:
                exonlist = dict_exonlist[ensg]
                exonstring = '|'.join(exonlist)
                e1,e2 = exons.split('-')
                pattern1 = re.compile(r'^{}\|{}\|'.format(e1,e2))  # ^E1.1|E2.3|
                pattern2 = re.compile(r'\|{}\|{}$'.format(e1,e2))  # |E1.1|E2.3$
                pattern3 = re.compile(r'\|{}\|{}\|'.format(e1,e2)) # |E1.1|E2.3|
                if re.search(pattern3,exonstring) or re.search(pattern2,exonstring) or re.search(pattern1,exonstring):   # as long as match one pattern, should be eliminated
                    continue
                else:
                    updated_valid.append(uid)
        print('reduce valid Neojunction from {} to {} because they are present in Ensembl db'.format(len(valid),len(updated_valid)))
        valid = updated_valid
    # consider add_control
    mean_add_list = []
    if add_control is not None:
        for i,(id_,control) in enumerate(add_control.items()):
            n_previous_valid = len(valid)
            if isinstance(control,pd.DataFrame):
                junction_to_mean = control.mean(axis=1).to_dict()
            elif isinstance(control,ad.AnnData):
                junction_to_mean = control.to_df().mean(axis=1).to_dict()
            else:
                raise Exception('control must be either in dataframe or anndata format')
            df['mean_add'] = df.index.map(junction_to_mean).fillna(value=0)
            df['diff_add'] = df['max'] - df['mean_add']
            df['cond_add'] = (df['mean_add'] < n_max) & (df['diff_add'] > t_min)
            mean_add_list.append(df['mean_add'])
            valid_add = df.loc[df['cond_add']].index.tolist()
            valid = list(set(valid).intersection(set(valid_add)))
            tmp = df.copy(); tmp.drop(columns=['mean','diff','cond'],inplace=True); tmp.rename(columns=lambda x:x+'_{}'.format(id_),inplace=True)
            df_to_write.append(tmp)
            print('reduce valid Neojunction from {} to {} because they are present in added control {}'.format(n_previous_valid,len(valid),id_))
    invalid = list(set(junction_count_matrix.index).difference(set(valid)))
    # now, consider each entry
    gtex_df = pd.concat([df['mean']]*junction_count_matrix.shape[1],axis=1)
    gtex_df.columns = junction_count_matrix.columns
    diff_df_gtex = junction_count_matrix - gtex_df
    cond_df = (gtex_df < n_max) & (diff_df_gtex > t_min)
    if add_control is not None:
        for mean_add in mean_add_list:
            add_df = pd.concat([mean_add]*junction_count_matrix.shape[1],axis=1)    
            add_df.columns = junction_count_matrix.columns
            diff_df_add = junction_count_matrix - add_df
            cond_df = cond_df & (add_df < n_max) & (diff_df_add > t_min)
    # write the df
    df_to_write = pd.concat(df_to_write,axis=1)
    df_to_write.to_csv(os.path.join(outdir,'NeoJunction_statistics_maxmin.txt'),sep='\t')
    return valid,invalid,cond_df


def crude_tumor_specificity(uid,count):    # for NeoJunction class, since we normally start from Jcmq with check_gtex=False, rarely being called.
    detail = ''
    if uid not in set(adata.obs_names):
        mean_value = 0
    else:
        mean_value = adata.obs.loc[uid,'mean']
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

def split_df_to_chunks(df,cores=None):
    df_index = np.arange(df.shape[0])
    if cores is None:
        cores = mp.cpu_count()
    sub_indices = np.array_split(df_index,cores)
    sub_dfs = [df.iloc[sub_index,:] for sub_index in sub_indices]
    return sub_dfs


def split_array_to_chunks(array,cores=None):
    if not isinstance(array,list):
        raise Exception('split_array_to_chunks function works for list, not ndarray')
    array_index = np.arange(len(array))
    if cores is None:
        cores = mp.cpu_count()
    sub_indices = np.array_split(array_index,cores)
    sub_arrays = []
    for sub_index in sub_indices:
        item_in_group = []
        for i in sub_index:
            item_in_group.append(array[i])
        sub_arrays.append(item_in_group)
    return sub_arrays


def add_tumor_specificity_frequency_table(df,method='mean',remove_quote=True,cores=None):
    '''
    add tumor specificty to each neoantigen-uid in the frequency table produced by SNAF T pipeline

    :param df: DataFrame, the frequency table produced by SNAF T pipeline
    :param method: string, either 'mean', or 'mle', or 'bayesian'
    :param remove quote: boolean, whether to remove the quotation or not, as one column in frequency table df is list, when loaded in memory using pandas, it will be added a quote, we can remove it
    :param cores: int, how many cpu cores to use for this computation, default None and use all the cpu the program detected

    :return new_df: a dataframe with one added column containing tumor specificity score

    Example::

        snaf.add_tumor_specificity_frequency_table(df,'mle',remove_quote=True)

    '''
    from ast import literal_eval
    import multiprocessing as mp
    if remove_quote:
        df['samples'] = [literal_eval(item) for item in df['samples']]
    if cores is None:
        cores = mp.cpu_count()

    if method != 'bayesian': 
        pool = mp.Pool(processes=cores)
        print('{} subprocesses have been spawned'.format(cores))

        all_unique_junctions = list(set([item.split(',')[1] for item in df.index]))
        sub_arrays = split_array_to_chunks(all_unique_junctions,cores=cores)
        r = [pool.apply_async(func=add_tumor_specificity_frequency_table_atomic_func,args=(sub_array,method,)) for sub_array in sub_arrays]  

        pool.close()
        pool.join()
        results = []
        for collect in r:
            result = collect.get()
            results.append(result)
        all_score_dict = {}
        for score_dict in results:
            all_score_dict.update(score_dict)
        col = []
        for item in df.index:
            col.append(all_score_dict[item.split(',')[1]])
        new_df = df.copy()
        new_df['tumor_specificity_{}'.format(method)] = col

    else:   # seems like bayesian doesn't work well with multiprocessing
        all_unique_junctions = list(set([item.split(',')[1] for item in df.index]))
        score_dict = {}
        for uid in tqdm(all_unique_junctions,total=len(all_unique_junctions)):
            score_dict[uid] = tumor_specificity(uid,'bayesian')
        for item in df.index:
            col.append(score_dict[item.split(',')[1]])
        new_df = df.copy()
        new_df['tumor_specificity_{}'.format(method)] = col
            

    return new_df

def add_tumor_specificity_frequency_table_atomic_func(sub_array,method):
    uid_list = sub_array
    score_dict = {uid:tumor_specificity(uid,method) for uid in tqdm(uid_list,total=len(uid_list))}
    return score_dict


    

def tumor_specificity(uid,method,return_df=False):
    try:
        info = adata[[uid],:]
    except:
        print('{} not detected in gtex, impute as zero'.format(uid))
        info = ad.AnnData(X=csr_matrix(np.full((1,adata.shape[1]),0)),obs=pd.DataFrame(data={'mean':[0]},index=[uid]),var=adata.var)  # weired , anndata 0.7.6 can not modify the X in place? anndata 0.7.2 can do that in scTriangulate
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
        scale_factor_dict = adata.var['total_count'].to_dict()
        df['value_cpm'] = df['value'].values / df.index.map(scale_factor_dict).values
        y = df['value_cpm'].values
        # mle_model = minimize(mle_func,np.array([0.2]),args=(y,),bounds=((0,1),),method='Nelder-Mead')
        mle_model = minimize_scalar(mle_func,bounds=(0,1),args=(y,),method='bounded')
        if mle_model.success:
            sigma = mle_model.x
        else:   
            sigma = 0
            print(uid,y, mle_model)   # debug purpose
        if return_df:
            return sigma,df
        else:
            return sigma
    elif method == 'bayesian':
        scale_factor_dict = adata.var['total_count'].to_dict()
        df['value_cpm'] = df['value'].values / df.index.map(scale_factor_dict).values
        y = df['value_cpm'].values
        x = []
        for tissue in adata.var['tissue'].unique():
            sub = adata[uid,adata.var['tissue']==tissue]
            total_count = sub.shape[1]
            c = np.count_nonzero(sub.X.toarray())
            scaled_c = round(c * (25/total_count),0)
            x.append(scaled_c)
        x = np.array(x)
        try:
            with pm.Model() as m:
                sigma = pm.Uniform('sigma',lower=0,upper=1)
                nc = pm.HalfNormal('nc',sigma=sigma,observed=y)
                nc_hat = pm.Deterministic('nc_hat',pm.math.sum(nc)/len(y))
                psi = pm.Beta('psi',alpha=2,beta=nc_hat*20)
                mu = pm.Gamma('mu',alpha=nc_hat*50,beta=1)
                c = pm.ZeroInflatedPoisson('c',psi,mu,observed=x)
                trace = pm.sample(draws=1000,step=pm.NUTS(target_accept=0.95),tune=1000,return_inferencedata=False,cores=1)
                '''
                the error of "Got error No model on context stack. trying to find log_likelihood in translation" maybe due to pymc build and how they launch multi-cores.
                remember, my build can only work when cores=1, which further indicate there might be an issue revolving around it.
                https://stackoverflow.com/questions/69888492/sampling-of-pymc3-in-python-gets-runtime-error-of-bootstrapping-phase
                '''

            df = az.summary(trace,round_to=2)

            '''
            az.summary(trace)

                    mean    sd  hdi_3%  hdi_97%  mcse_mean  mcse_sd  ess_bulk  ess_tail  r_hat
            sigma    0.47  0.01    0.46     0.48       0.00     0.00    182.23     98.84   1.02
            nc_hat   0.22  0.00    0.22     0.22       0.00     0.00    200.00    200.00    NaN
            mu      22.97  0.52   21.87    23.87       0.04     0.03    196.69     94.84   1.00

            az.plot_posterior(trace,var_names=['sigma','nc_hat','mu'])
            az.plot_forest(trace,,var_names=['sigma','nc_hat','mu'])

            gv = pm.model_to_graphviz(m)
            gv.format = 'pdf'
            gv.render(filename='model_graph');sys.exit('stop')
            # to run the above, you need to module load graphviz so that dot is exposed to the program

            '''
            sigma = df.iloc[0]['mean']

        except:
            sigma = None
            print(uid,x,y)

        if return_df:
            return sigma,df
        else:
            return sigma
            
















