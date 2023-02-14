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
from kneed import KneeLocator


# functions
def compute_y(adata,uids):
    info = adata[uids,:]
    y = info.X.toarray() / adata.var['total_count'].values.reshape(1,-1)
    return y

def compute_x(adata,uids):
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

def thresholding_kneedle(cpm,plot=False):
    x = cpm[cpm > 0]  # all non-zero values
    actual_x = np.arange(len(x))
    actual_y = np.sort(x)
    try:
        kneedle = KneeLocator(actual_x,actual_y,S=1,curve='convex',direction='increasing',interp_method='polynomial')
    except:
        
    knee_index = round(kneedle.knee,0)
    if plot:
        kneedle.plot_knee()
        plt.savefig('kneedle_data.pdf',bbox_inches='tight')
        plt.close()
        kneedle.plot_knee_normalized()
        plt.savefig('kneedle_norm.pdf',bbox_inches='tight')
        plt.close()
    return actual_y[knee_index]

def threshold(cpm,method,**kwargs):
    if method == 'kneedle':
        th = thresholding_kneedle(cpm,**kwargs)
    elif method == 'otsu':
        th = thresholding_otsu(cpm,**kwargs)
    cpm = np.where(cpm>th,cpm,0)
    return cpm

def thresholding_otsu(cpm,step=0.05,dampen_factor=20):
    x = cpm[cpm > 0]  # all non-zero values
    criteria = []
    ths = np.arange(0,x.max(),step)
    for th in ths:
        thresholded_x = np.where(x>=th,1,0)
        w1 = np.count_nonzero(thresholded_x)/len(x)
        w0 = 1 - w1
        if w1 == 0 or w0 == 0:
            value = np.inf
        else:
            x1 = x[thresholded_x==1]
            x0 = x[thresholded_x==0]
            var1 = x1.var()
            var0 = x0.var()
            value = w0 * var0 + w1 * var1 / dampen_factor
        criteria.append(value)
    best_th = ths[np.argmin(criteria)]
    return best_th


adata = ad.read_h5ad('coding.h5ad')
uids = adata.obs_names.tolist()

# from gtex_viewer import *
# gtex_viewer_configuration(adata)
# df = gtex_visual_norm_count_combined('ENSG00000198681',xlim=None,ylim=(0,10),save_df=False)
Y = compute_y(adata,uids)
thresholded_Y = np.empty_like(Y,dtype=np.float64)
for i in tqdm(range(Y.shape[0]),total=Y.shape[0]):
    thresholded_Y[i,:] = threshold(Y[i,:],'kneedle')
sys.exit('stop')


X = compute_scaled_x(adata,uids)



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



