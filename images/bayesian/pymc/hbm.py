import anndata as ad
import numpy as np
import pandas as pd
import os,sys
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle

import arviz as az
import pymc as pm
import pytensor.tensor as at



def compute_scaled_x(adata,uid):
    x = []
    for tissue in adata.var['tissue'].unique():
        sub = adata[uid,adata.var['tissue']==tissue]
        total_count = sub.shape[1]
        if total_count >= 10:
            c = np.count_nonzero(np.where(sub.X.toarray()<1,0,sub.X.toarray()))
            scaled_c = round(c * (25/total_count),0)
            x.append(scaled_c)
    x = np.array(x)
    return x

def compute_y(adata,uid):
    info = adata[[uid],:]
    scale_factor_dict = adata.var['total_count'].to_dict()
    df = pd.DataFrame(data={'value':info.X.toarray().squeeze(),'tissue':info.var['tissue'].values},index=info.var_names)
    df['value_cpm'] = df['value'].values / df.index.map(scale_factor_dict).values
    y = df['value_cpm'].values 
    return y

def diagnose(final_path,output_name='diagnosis.pdf'):
    df = pd.read_csv(final_path,sep='\t',index_col=0)
    fig,ax = plt.subplots()
    im = ax.scatter(df['mean_x'],df['mean_y'],c=df['sigma'])
    plt.colorbar(im)
    ax.set_ylabel('average_normalized_counts')
    ax.set_xlabel('average_n_present_samples_per_tissue')
    plt.savefig(output_name,bbox_inches='tight')
    plt.close()

def prior_check(prior_samples,var,observed):
    fig,ax = plt.subplots()
    az.plot_dist(
        observed,
        kind="hist",
        color="orange",
        hist_kwargs=dict(alpha=0.6,bins=100),
        label="observed",
        ax = ax
    )
    az.plot_dist(
        prior_samples.prior_predictive[var],
        kind="hist",
        color = 'blue',
        hist_kwargs=dict(alpha=0.6,bins=100),
        label="simulated",
        ax = ax
    )
    ax.tick_params(axis='x',rotation=45,labelsize=0.5)
    plt.savefig('prior_check_{}.pdf'.format(var),bbox_inches='tight')
    plt.close()

def posterior_check(posterior_samples,var,observed):
    fig,ax = plt.subplots()
    az.plot_dist(
        observed,
        kind="hist",
        color="orange",
        hist_kwargs=dict(alpha=0.6,bins=100),
        label="observed",
        ax = ax
    )
    az.plot_dist(
        posterior_samples.posterior_predictive[var],
        kind="hist",
        color = 'blue',
        hist_kwargs=dict(alpha=0.6,bins=100),
        label="simulated",
        ax = ax
    )
    ax.tick_params(axis='x',rotation=45,labelsize=0.5)
    plt.savefig('posterior_check_{}.pdf'.format(var),bbox_inches='tight')
    plt.close() 

def infer_parameters_vectorize(uids,solver='mcmc'):
    Y = np.array([compute_y(adata,uid) for uid in uids])
    X = np.array([compute_scaled_x(adata,uid) for uid in uids])
    n = len(uids)
    s = Y.shape[1]
    t = X.shape[1]
    with pm.Model() as m:
        sigma = pm.Uniform('sigma',lower=0,upper=1,shape=n)
        nc = pm.HalfNormal('nc',sigma=sigma,observed=Y.T)  # since sigma is of shape n, so each draw will be a row of length n (support dimention), and we draw s times, then the resultant matrix, we transpose
        psi = pm.Beta('psi',alpha=2,beta=sigma*2)
        mu = pm.Gamma('mu',alpha=sigma*25,beta=1)
        c = pm.ZeroInflatedPoisson('c',psi,mu,observed=X.T)
        # gv = pm.model_to_graphviz(m)
        # gv.format = 'pdf'
        # gv.render(filename='model_graph')
    with m:
        if solver == 'mcmc':
            trace = pm.sample(draws=1000,step=pm.NUTS(),tune=1000,cores=1,progressbar=True)
        elif solver == 'vi':
            mean_field = pm.fit(method='advi',progressbar=False)
            trace = mean_field.sample(1000)
    with open('pickle_vi_trace.p','wb') as f:
        pickle.dump(trace,f)

def infer_parameters(uid):
    y = compute_y(adata,uid)
    x = compute_scaled_x(adata,uid)
    with pm.Model() as m:
        sigma = pm.Uniform('sigma',lower=0,upper=1)
        nc = pm.HalfNormal('nc',sigma=sigma,observed=y)
        psi = pm.Beta('psi',alpha=2,beta=sigma*2)  # beta(2,2) is good for modeling intermediate proportion, beta(0.5,0.5) is good for modeling bimodel
        mu = pm.Gamma('mu',alpha=sigma*25,beta=1)  # gamme(alpha,1) is good to control the mean, as mean=alpha/beta
        c = pm.ZeroInflatedPoisson('c',psi,mu,observed=x)
    # gv = pm.model_to_graphviz(m)
    # gv.format = 'pdf'
    # gv.render(filename='model_graph')
    # with m:
    #     prior_samples = pm.sample_prior_predictive(1000)
    # prior_check(prior_samples,'nc',y)
    # prior_check(prior_samples,'c',x)
    # with m:
    #     trace = pm.sample(draws=1000,step=pm.NUTS(),tune=1000,cores=1,progressbar=False)
    # df = az.summary(trace,round_to=4)
    # result = df['mean'].tolist() + [y.mean(),np.array(x).mean()]
    with m:
        # below is from pymc3 tutorial: https://docs.pymc.io/en/v3/pymc-examples/examples/variational_inference/variational_api_quickstart.html
        mean_field = pm.fit(method='advi',progressbar=False)
    posterior_samples = mean_field.sample(1000).posterior
    result = [posterior_samples['sigma'].values.mean(),posterior_samples['psi'].values.mean(),posterior_samples['mu'].values.mean()] + [y.mean(),np.array(x).mean()]
    # az.plot_trace(trace,var_names=['sigma','mu','psi'])
    # plt.savefig('trace.pdf',bbox_inches='tight');plt.close()
    # az.plot_energy(trace)
    # plt.savefig('energy.pdf',bbox_inches='tight');plt.close()
    # az.plot_forest(trace,var_names=['sigma','mu','psi'],combined=True,hdi_prob=0.95,r_hat=True)
    # plt.savefig('forest.pdf',bbox_inches='tight');plt.close()
    # with m:
    #     extended_trace = pm.sample_posterior_predictive(trace)
    # posterior_check(extended_trace,'c',x)
    # posterior_check(extended_trace,'nc',y)
    return result


# load the count
# adata = ad.read_h5ad('combined_normal_count.h5ad')
# adata.var['tissue'].value_counts().to_csv('tissue_count.txt',sep='\t')
# adata = adata[np.random.choice(np.arange(adata.shape[0]),size=100,replace=False),:]
# adata.write('sampled_100.h5ad')
# adata.to_df().to_csv('sampled_100.txt',sep='\t')
adata = ad.read_h5ad('sampled_100.h5ad')

# # visualize
# uid = 'ENSG00000115459:I5.1-E6.1'
# from gtex_viewer import *
# gtex_viewer_configuration(adata)
# gtex_visual_per_tissue_count(uid)
# gtex_visual_norm_count_combined(uid)
# gtex_visual_subplots(uid)

# # determine the coefficient between n_hat and psi/mu, but seems to be indirect to model nc_hat and psi/mu
# X = np.array([compute_scaled_x(adata,uid) for uid in adata.obs_names])
# Y = np.array([compute_y(adata,uid) for uid in adata.obs_names]).mean(axis=1)
# psi = np.count_nonzero(X,axis=1) / X.shape[1]
# mu = np.quantile(X,0.5,axis=1)
# nc_hat = Y
# df = adata.to_df()
# df['psi'] = psi; df['mu'] = mu; df['nc_hat'] = nc_hat
# df.to_csv('empirical.txt',sep='\t')
# train_X = np.column_stack([np.ones(len(nc_hat)),nc_hat])
# # linear regression
# result = np.linalg.lstsq(X,np.column_stack([psi,mu]),rcond=None)
# # GLM
# import statsmodels.api as sm
# mod = sm.GLM(endog=psi,exog=train_X,family=sm.families.Binomial())
# mod_result = mod.fit()
# mod = sm.GLM(endog=mu,exog=train_X,family=sm.families.Poisson())
# mod_result = mod.fit()
# mod_result.summary()
# fig,ax = plt.subplots()
# ax.scatter(nc_hat,mu)
# plt.savefig('mu.pdf',bbox_inches='tight')
# plt.close()

# change to model
# uid = 'ENSG00000115459:I5.1-E6.1'
# infer_parameters(uid)

# # mode
# count = pd.read_csv('sampled_100.txt',sep='\t',index_col=0)
# data = []
# for uid in tqdm(adata.obs_names,total=adata.shape[0]):
#     values = infer_parameters(uid)
#     data.append((uid,*values))
# result = pd.DataFrame.from_records(data,columns=['uid','sigma','psi','mu','mean_y','mean_x']).set_index('uid')
# final = pd.concat([count,result],axis=1)
# final.to_csv('final.txt',sep='\t')

# # diagnose
# diagnose('final.txt','diagnosis.pdf')

# vectorize
# infer_parameters_vectorize(uids=adata.obs_names.tolist(),solver='vi')
count = pd.read_csv('sampled_100.txt',sep='\t',index_col=0)
with open('pickle_vi_trace.p','rb') as f:
    trace = pickle.load(f)
df = az.summary(trace,round_to=4)
values_list = []
for param in ['sigma','psi','mu']:
    tmp = df.loc[df.index.to_series().str.contains(pat=param),['mean','r_hat']].set_index(count.index).rename(columns=lambda x:x+'_{}'.format(param))
    values_list.append(tmp)
final = pd.concat([count,*values_list],axis=1)
final.to_csv('final_vectorize_vi.txt',sep='\t')












