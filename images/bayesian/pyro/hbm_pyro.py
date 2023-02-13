#!/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ts/bayesian/pytorch_pyro_mamba_env/bin/python3.7

import anndata as ad  # need to install from -c bioconda, not -c conda-forge, but still fail (install 0.6.2 instead), finally use pip to solve (0.8.0)
import numpy as np
import pandas as pd
import os,sys
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import pickle
import torch
import pyro
import pickle
import pyro.poutine as poutine
import pyro.distributions as dist
import pyro.distributions.constraints as constraints
from pyro.infer import SVI,Trace_ELBO
from pyro.optim import Adam,ClippedAdam

# functions
def compute_scaled_x(adata,uids):
    total_tissue = adata.var['tissue'].unique()
    valid_tissue = [tissue for tissue in total_tissue if adata[:,adata.var['tissue']==tissue].shape[1] >= 10]
    x = np.zeros((len(uids),len(valid_tissue)))
    for i,tissue in enumerate(valid_tissue):
        sub = adata[uids,adata.var['tissue']==tissue]
        total_count = sub.shape[1]
        c = np.count_nonzero(np.where(sub.X.toarray()<1,0,sub.X.toarray()),axis=1)
        scaled_c = np.round(c * (25/total_count),0)
        x[:,i] = scaled_c
    return x

def compute_y(adata,uids):
    info = adata[uids,:]
    y = info.X.toarray() / adata.var['total_count'].values.reshape(1,-1)
    return y

def diagnose(final_path,output_name='diagnosis.pdf'):
    df = pd.read_csv(final_path,sep='\t',index_col=0)
    fig,ax = plt.subplots()
    im = ax.scatter(df['X_mean'],df['Y_mean'],c=df['mean_sigma'],s=0.5**2,cmap='viridis')
    plt.colorbar(im)
    ax.set_ylabel('average_normalized_counts')
    ax.set_xlabel('average_n_present_samples_per_tissue')
    ax.set_ylim([-1,7])
    plt.savefig(output_name,bbox_inches='tight')
    plt.close()

def diagnose_plotly(final_path,output_name='diagnosis.html'):
    import plotly.graph_objects as go
    df = pd.read_csv(final_path,sep='\t',index_col=0)
    df['symbol'] = ensemblgene_to_symbol(df.index.tolist(),'human')
    df.loc[:,['X_mean','Y_mean','mean_sigma','symbol']].to_csv('result.txt',sep='\t')
    node_x = []
    node_y = []
    node_text = []
    node_color = []
    for uid,X_mean,Y_mean,mean_sigma in zip(df.index,df['X_mean'],df['Y_mean'],df['mean_sigma']):
        if Y_mean < 7:
            node_x.append(X_mean)
            node_y.append(Y_mean)
            node_text.append('{};X:{};Y:{}'.format(uid,X_mean,Y_mean))
            node_color.append(mean_sigma)
    node_trace = go.Scatter(x=node_x,y=node_y,mode='markers',marker={'color':node_color,'colorscale':'Viridis','showscale':True,'size':1},text=node_text,hoverinfo='text')
    fig_layout = go.Layout(showlegend=False,title='diagnose',xaxis=dict(title_text='average_n_present_samples_per_tissue'),yaxis=dict(title_text='average_normalized_counts'))
    fig = go.Figure(data=[node_trace],layout=fig_layout)
    fig.write_html(output_name,include_plotlyjs='cdn')

# test
adata = ad.read_h5ad('coding.h5ad')
uids = adata.obs_names.tolist()
# X = compute_scaled_x(adata,uids)
# Y = compute_y(adata,uids)
# with open('X.p','wb') as f:
#     pickle.dump(X,f)
# with open('Y.p','wb') as f:
#     pickle.dump(Y,f)

with open('X.p','rb') as f:
    X = pickle.load(f)
with open('Y.p','rb') as f:
    Y = pickle.load(f)
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device)
X = torch.tensor(X.T,device=device)
Y = torch.tensor(Y.T,device=device)
n = X.shape[1]
s = Y.shape[0]
t = X.shape[0]


adam = Adam({'lr': 0.002,'betas':(0.95,0.999)}) 
clipped_adam = ClippedAdam({'betas':(0.95,0.999)})
elbo = Trace_ELBO()

def model_mle(X,Y):
    sigma = pyro.param('sigma',lambda:torch.tensor(np.full(n,0.5),device=device),constraint=constraints.interval(0.05,1))
    with pyro.plate('data_Y'):
        nc = pyro.sample('nc',dist.HalfNormal(sigma).expand([s,n]).to_event(1),obs=Y)
    psi = pyro.sample('psi',dist.Beta(sigma*20.,torch.tensor(10.,device=device)).to_event(1))
    mu = pyro.sample('mu',dist.Gamma(sigma*25.,torch.tensor(1.,device=device)).to_event(1))
    with pyro.plate('data_X'):
        c = pyro.sample('c',dist.ZeroInflatedPoisson(rate=mu,gate=1.-psi).expand([t,n]).to_event(1),obs=X)

def guide_mle(X,Y):
    pass

# train
n_steps = 2000
pyro.clear_param_store()
svi = SVI(model_mle, guide_mle, adam, loss=Trace_ELBO())
losses = []
for step in tqdm(range(n_steps),total=n_steps):  
    loss = svi.step(X,Y)
    losses.append(loss)
    print("Elbo loss step {}: {}".format(step,loss))
plt.figure(figsize=(5, 2))
plt.plot(losses)
plt.xlabel("SVI step")
plt.ylabel("ELBO loss")
plt.savefig('elbo_loss.pdf',bbox_inches='tight')
plt.close()

sigma = pyro.param('sigma').data.cpu().numpy()
df = pd.Series(index=uids,data=sigma,name='mean_sigma').to_frame()
Y_mean = Y.mean(axis=0).data.cpu().numpy()
X_mean = X.mean(axis=0).data.cpu().numpy()
df['Y_mean'] = Y_mean
df['X_mean'] = X_mean
df.to_csv('mle_results.txt',sep='\t')
diagnose('mle_results.txt','pyro_mle_diagnosis.pdf')
sys.exit('stop')



def model(X,Y):
    sigma = pyro.sample('sigma',dist.Uniform(torch.tensor(0.,device=device),torch.tensor(1.,device=device)).expand([n]).to_event(1))
    with pyro.plate('data_Y'):
        nc = pyro.sample('nc',dist.HalfNormal(sigma).expand([s,n]).to_event(1),obs=Y)
    psi = pyro.sample('psi',dist.Beta(sigma*20.,torch.tensor(10.,device=device)).to_event(1))
    mu = pyro.sample('mu',dist.Gamma(sigma*25.,torch.tensor(1.,device=device)).to_event(1))
    with pyro.plate('data_X'):
        c = pyro.sample('c',dist.ZeroInflatedPoisson(rate=mu,gate=1.-psi).expand([t,n]).to_event(1),obs=X)

# pyro.render_model(model, model_args=(X,Y), render_distributions=True, filename='bayesian_pyro.pdf')
trace = poutine.trace(model).get_trace(X,Y)
trace.compute_log_prob()  
print(trace.format_shapes())





def guide(X,Y):
    sigma_alpha = pyro.param('sigma_alpha',lambda: torch.tensor(np.full(shape=n,fill_value=2.),device=device),constraint=constraints.positive)
    sigma = pyro.sample('sigma',dist.Beta(sigma_alpha,torch.tensor([2],device=device)).expand([n]).to_event(1))

    psi_alpha = pyro.param('psi_alpha',lambda: torch.tensor(np.full(shape=n,fill_value=10.),device=device),constraint=constraints.positive)
    psi_beta = pyro.param('psi_beta',lambda: torch.tensor(np.full(shape=n,fill_value=10.),device=device),constraint=constraints.positive)
    psi = pyro.sample('psi',dist.Beta(psi_alpha,psi_beta).expand([n]).to_event(1))

    mu_alpha = pyro.param('mu_mean',lambda: torch.tensor(np.full(shape=n,fill_value=12.5),device=device),constraint=constraints.positive)
    mu_beta = pyro.param('mu_scale',lambda: torch.ones(n,device=device),constraint=constraints.positive)
    mu = pyro.sample('mu',dist.Gamma(mu_alpha,mu_beta).expand([n]).to_event(1))
    return {'sigma':sigma,'psi':psi,'mu':mu}

# guide = pyro.infer.autoguide.AutoNormal(model)
# pyro.render_model(guide, model_args=(X,Y), render_distributions=True, render_params=True, filename='guide.pdf')
trace = poutine.trace(guide).get_trace(X,Y)
trace.compute_log_prob()  
print(trace.format_shapes())

pyro.clear_param_store()
svi = pyro.infer.SVI(model, guide, adam, elbo)

losses = []
for step in range(1000):  
    loss = svi.step(X,Y)
    losses.append(loss)
    print("Elbo loss: {}".format(loss))
plt.figure(figsize=(5, 2))
plt.plot(losses)
plt.xlabel("SVI step")
plt.ylabel("ELBO loss")
plt.savefig('elbo_loss.pdf',bbox_inches='tight')
plt.close()
sys.exit('stop')

for name, value in pyro.get_param_store().items():
    print(name, pyro.param(name).data.cpu().numpy())
    if name == 'sigma':
        pd.Series(data=pyro.param(name).data.cpu().numpy())

predictive = pyro.infer.Predictive(model, guide=auto_guide, num_samples=10000)
svi_samples = predictive(X,Y) # same dict as samples above, but here they have obs, obs is generated by first sample from guide, then run model using those values
svi_nc = svi_samples['nc'] # (800,170)
svi_c = svi_samples['c'] # (800,170)














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

def infer_parameters_vectorize_pyro(uids):
    X = compute_scaled_x(adata,uids)
    Y = compute_y(adata,uids)
    n = len(uids)
    s = Y.shape[1]
    t = X.shape[1]

    import pyro.distributions as dist
    import pyro.distributions.constraints as constraints
    X = torch.tensor(X)
    Y = torch.tensor(Y)
    def model(X,Y):
        sigma = pyro.sample('sigma',dist.Uniform(0,1).expand([n]))
        with pyro.plate('data_Y'):
            nc = pyro.sample('nc',dist.HalfNormal(sigma).expand([s,n]),obs=Y.T)
        psi = pyro.sample('psi',dist.Beta(2,sigma*2))
        mu = pyro.sample('mu',dist.Gamma(sigma*25,1))
        with pyro.plate('data_X'):
            c = pyro.sample('c',dist.ZeroInflatedPoisson(rate=mu,gate=1-psi).expand([t,n]),obs=X.T)
        return nc,c
    pyro.render_model(model, model_args=(X,Y), render_distributions=True, filename='bayesian_pyro.pdf')
    sys.exit('stop')

def infer_parameters_vectorize(uids,solver='vi'):
    X = compute_scaled_x(adata,uids)
    Y = compute_y(adata,uids)
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
            mean_field = pm.fit(method='advi',progressbar=True)
            trace = mean_field.sample(1000)
    with open('pickle_trace.p','wb') as f:
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
        mean_field = pm.fit(method='advi',progressbar=True)
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

def ensemblgene_to_symbol(query,species):
    # assume query is a list, will also return a list
    import mygene
    mg = mygene.MyGeneInfo()
    out = mg.querymany(query,scopes='ensemblgene',fileds='symbol',species=species,returnall=True,as_dataframe=True,df_index=True)

    df = out['out']
    df_unique = df.loc[~df.index.duplicated(),:]
    df_unique['symbol'].fillna('unknown_gene',inplace=True)
    mapping = df_unique['symbol'].to_dict()

    result = []
    for item in query:
        result.append(mapping[item])

    return result

# # load the count
# adata = ad.read_h5ad('combined_gene_count.h5ad')
# genes = ensemblgene_to_symbol(query=adata.obs_names.tolist(),species='human')
# adata.obs['symbol'] = genes
# adata.obs.to_csv('mean_and_symbol_inspection.txt',sep='\t')

# # summarize LINC, LOC and -AS genes versus coding gene
# df = pd.read_csv('mean_and_symbol_inspection.txt',sep='\t',index_col=0)
# linc = df.loc[df['symbol'].str.startswith('LINC'),:]
# loc = df.loc[df['symbol'].str.startswith('LOC'),:]
# antisense = df.loc[df['symbol'].str.contains('-AS'),:]
# unknown = df.loc[df['symbol']=='unknown_gene',:]
# cond = ~((df['symbol'].str.startswith('LINC')) | (df['symbol'].str.startswith('LOC')) | (df['symbol'].str.contains('-AS')) | (df['symbol']=='unknown_gene'))
# coding = df.loc[cond,:]
# fig,ax = plt.subplots()
# sns.boxplot(data=[linc['mean'],loc['mean'],antisense['mean'],unknown['mean'],coding['mean']],ax=ax)
# ax.set_xticklabels(['LINC','LOC','Antisense','unknown','coding gene'])
# ax.set_ylabel('raw count')
# plt.savefig('summary_category.pdf',bbox_inches='tight')
# plt.close()

# adata = adata[coding.index,:]
# adata.write('coding.h5ad')  # 24290 × 3644

# adata.var['tissue'].value_counts().to_csv('tissue_count.txt',sep='\t')
# adata = adata[np.random.choice(np.arange(adata.shape[0]),size=20000,replace=False),:]
# adata.write('sampled_20000.h5ad')
# adata.to_df().to_csv('sampled_20000.txt',sep='\t')
# adata = ad.read_h5ad('sampled_20000.h5ad')

adata = ad.read_h5ad('coding.h5ad')
# adata.to_df().to_csv('coding.txt',sep='\t')


# # visualize
# from gtex_viewer import *
# gtex_viewer_configuration(adata)
# # selected = adata[np.random.choice(np.arange(adata.shape[0]),size=100,replace=False),:].obs_names.tolist()
# def visualize(uid):
#     outdir = 'random_100_dist/{}'.format(uid)
#     gtex_visual_per_tissue_count(uid,outdir=outdir)
#     gtex_visual_norm_count_combined(uid,outdir=outdir)
#     gtex_visual_subplots(uid,outdir=outdir,top=20)
# # for uid in tqdm(selected):
# #     visualize(uid)
# uid = 'ENSG00000185664'
# visualize(uid)


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

# vectorize
uids = adata.obs_names.tolist()
infer_parameters_vectorize_pyro(uids=uids)
sys.exit('stop')

count = pd.read_csv('coding.txt',sep='\t',index_col=0)
with open('pickle_trace.p','rb') as f:
    trace = pickle.load(f)
df = az.summary(trace,round_to=4)
values_list = []
for param in ['sigma','psi','mu']:
    tmp = df.loc[df.index.to_series().str.contains(pat=param),['mean','r_hat']].set_index(count.index).rename(columns=lambda x:x+'_{}'.format(param))
    values_list.append(tmp)
final = pd.concat([count,*values_list],axis=1)
uids = adata.obs_names.tolist()
Y_mean = compute_y(adata,uids).mean(axis=1)
X_mean = compute_scaled_x(adata,uids).mean(axis=1)
final['Y_mean'] = Y_mean
final['X_mean'] = X_mean
final.to_csv('final_vectorize_vi.txt',sep='\t')

# diagnose
diagnose_plotly('final_vectorize_vi.txt','diagnosis.html')












