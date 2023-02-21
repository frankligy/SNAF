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
from scipy.sparse import csr_matrix
from pyro.poutine import scale
from sklearn.mixture import GaussianMixture
from scipy.stats import pearsonr, spearmanr
import numpy.ma as ma
from sklearn.metrics import precision_recall_curve,auc,roc_curve,accuracy_score

# functions
def compute_y(adata,uids):
    info = adata[uids,:]
    y = info.X.toarray() / adata.var['total_count'].values.reshape(1,-1)
    return y

def compute_x(adata,uids):
    total_tissue = adata.var['tissue'].unique()
    valid_tissue = [tissue for tissue in total_tissue if adata[:,adata.var['tissue']==tissue].shape[1] >= 10]
    pd.Series(data=valid_tissue).to_csv('valid_tissue.txt',sep='\t')
    x = np.zeros((len(uids),len(valid_tissue)))
    for i,tissue in enumerate(valid_tissue):
        sub = adata[uids,adata.var['tissue']==tissue]
        total_count = sub.shape[1]
        c = np.count_nonzero(np.where(sub.X.toarray()<1,0,sub.X.toarray()),axis=1) / total_count
        x[:,i] = c
    return x

def get_thresholded_adata(adata,cond_Y):
    adata = adata.copy()
    adata.X = csr_matrix(np.where(cond_Y,adata.X.toarray(),0))
    return adata

def weighting(adata,dic,weights):
    total_tissue = adata.var['tissue'].unique()
    valid_tissue = [tissue for tissue in total_tissue if adata[:,adata.var['tissue']==tissue].shape[1] >= 10]
    for t,w in dic.items():
        i = valid_tissue.index(t)
        weights[i] = w
    return weights

def diagnose(final_path,ylim=(-1,200),output_name='diagnosis.pdf'):
    df = pd.read_csv(final_path,sep='\t',index_col=0)
    fig,ax = plt.subplots()
    im = ax.scatter(df['X_mean'],df['Y_mean'],c=df['mean_sigma'],s=0.5**2,cmap='viridis')
    plt.colorbar(im)
    ax.set_ylabel('average_normalized_counts')
    ax.set_xlabel('average_n_present_samples_per_tissue')
    ax.set_ylim(ylim)
    plt.savefig(output_name,bbox_inches='tight')
    plt.close()

def thresholding_kneedle(cpm,plot=False,S=1,interp_method='polynomial'):
    x = cpm[cpm > 0]  # all non-zero values
    if len(x) <= 2:
        return 0
    else:
        actual_x = np.arange(len(x))
        actual_y = np.sort(x)
        kneedle = KneeLocator(actual_x,actual_y,S=S,curve='convex',direction='increasing',interp_method=interp_method)
        knee = kneedle.knee
        if knee is not None:
            knee_index = round(knee,0)
        else:
            knee_index = 0
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
    elif method == 'gmm':
        th = thresholding_gmm(cpm,**kwargs)
    elif method == 'hardcode':
        th = thresholding_hardcode(cpm,**kwargs)
    cpm = np.where(cpm>th,cpm,0)
    cond = cpm>th
    return cpm, cond, th

def thresholding_hardcode(cpm,v):
    return v

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


def thresholding_gmm(cpm):
    x = cpm[cpm > 0]  # all non-zero values
    gm = GaussianMixture(n_components=2).fit(cpm.reshape(-1,1))
    means = gm.means_
    bg_index = np.argmin(means.mean(axis=1))
    best_th = means[bg_index,0]
    return best_th

def compute_concordance(uids, X,lookup,external,valid_tissue,method):
    ''' 
    X is n_gene * n_valid_tissue
    '''
    concordance = {}
    external_dic = {gene: sub_df for gene, sub_df in external.groupby(by='Gene')}   # {ensg:sub_df}
    external_ordered_tissue = {t:i for i, t in external_dic['ENSG00000121410']['Tissue'].reset_index(drop=True).to_dict().items()}  # {thyroid gland:0}
    external_dic = {gene: sub_df['nTPM'].values for gene, sub_df in external_dic.items()}
    lookup_dic = lookup.to_dict()  # {internal:external}
    valid_tissue = {t:i for i, t in valid_tissue.to_dict().items()}   # {thyroid: 4}
    avail_indices = []
    avail_external_indices = []
    for inter, exter in lookup_dic.items():
        avail_indices.append(valid_tissue[inter])
        avail_external_indices.append(external_ordered_tissue[exter])
    for i in tqdm(np.arange(X.shape[0]),total=X.shape[0]):
        gene = uids[i]
        v1 = X[i,avail_indices]
        try:
            v2 = external_dic[gene][avail_external_indices]
        except:
            continue
        if np.count_nonzero(v2) < len(v2) * 0.5 and np.count_nonzero(v1) > 0 and np.count_nonzero(v2) > 0:
            if method == 'pearson':
                value = pearsonr(v1,v2)[0]
            elif method == 'spearman':
                value = spearmanr(v1,v2)[0]
            elif method == 'AUPR':
                v2 = np.where(v2>0,1,0)
                precision,recall,_ = precision_recall_curve(v2,v1,pos_label=1)
                value = auc(recall,precision)
            concordance[gene] = value
    return concordance





'''main program starts'''

adata = ad.read_h5ad('coding.h5ad')
uids = adata.obs_names.tolist()

'''
I want to fix the threshold issue
testing targets:

ENSG00000198681  MAGEA1
ENSG00000221867  MAGEA3
ENSG00000156738 MS4A1, CD20
ENSG00000185274 GALNT17, threshould is about 30
ENSG00000141506, PIK3R5, threshold is about 20
ENSG00000150991 a highly expressed one
'''

# uid = 'ENSG00000156738'
# from gtex_viewer import *
# gtex_viewer_configuration(adata)
# df = gtex_visual_norm_count_combined(uid,xlim=None,ylim=None,save_df=True)
# gtex_visual_per_tissue_count(uid)
# gtex_visual_subplots(uid)

# lookup = pd.read_csv('tissue_lookup.txt',sep='\t',index_col=0).squeeze()
# external = pd.read_csv('rna_tissue_consensus.tsv',sep='\t')
# valid_tissue = pd.read_csv('valid_tissue.txt',sep='\t',index_col=0).squeeze()

# data = {}
# for v in np.arange(0,5.1,0.1):
#     print(v)
#     Y = compute_y(adata,uids)
#     thresholded_Y = np.empty_like(Y,dtype=np.float32)
#     cond_Y = np.empty_like(Y,dtype=bool)
#     ths = []
#     for i in tqdm(range(Y.shape[0]),total=Y.shape[0]):
#         thresholded_Y[i,:], cond_Y[i,:], th = threshold(Y[i,:],'hardcode',v=v)
#         ths.append(th)
#     # pd.Series(index=uids,data=ths,name='threshold').to_csv('threshold.txt',sep='\t')
#     new_adata = get_thresholded_adata(adata,cond_Y)
#     X = compute_x(new_adata,uids)
#     concordance = compute_concordance(uids,X,lookup,external,valid_tissue,'spearman')
#     print(ma.masked_invalid(list(concordance.values())).mean())
#     data[v] = concordance
# with open('concordance.p','wb') as f:
#     pickle.dump(data,f)

# with open('concordance.p','rb') as f:
#     data = pickle.load(f)
# x = list(data.keys())
# y = [ma.masked_invalid(list(v.values())).mean() for v in data.values()]
# # y_95 = [np.nanpercentile(list(v.values()),95) for v in data.values()]
# # y_05 = [np.nanpercentile(list(v.values()),5) for v in data.values()]
# fig,ax = plt.subplots()
# ax.plot(x,y,marker='o',linewidth=1,linestyle='-')
# # ax.fill_between(x,y_05,y_95,alpha=0.2)
# ax.set_xticks(x)
# ax.set_xticklabels(x,fontsize=0.5,rotation=90)
# ax.set_xlabel('cutoff')
# ax.set_ylabel('mean_spearmanr')
# plt.savefig('cutoff.pdf',bbox_inches='tight')
# plt.close()



'''
Seems 0.8 is a decent cutoff
'''
Y = compute_y(adata,uids)
thresholded_Y = np.empty_like(Y,dtype=np.float32)
cond_Y = np.empty_like(Y,dtype=bool)
ths = []
for i in tqdm(range(Y.shape[0]),total=Y.shape[0]):
    thresholded_Y[i,:], cond_Y[i,:], th = threshold(Y[i,:],'hardcode',v=0.8)
    ths.append(th)
# pd.Series(index=uids,data=ths,name='threshold').to_csv('threshold.txt',sep='\t')
new_adata = get_thresholded_adata(adata,cond_Y)
X = compute_x(new_adata,uids)

with open('raw_X.p','wb') as f:
    pickle.dump(X, f)

with open('raw_X.p','rb') as f:
    X = pickle.load(f)
with open('Y.p','rb') as f:
    Y = pickle.load(f)


# index = uids.index(uid)
# X = X[[index],:]
# Y = Y[[index],:]
# print(X.shape,Y.shape)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(device)
Y = np.where(Y==0,1e-5,Y)
X = torch.tensor(X.T,device=device)
Y = torch.tensor(Y.T,device=device)
n = X.shape[1]
s = Y.shape[0]
t = X.shape[0]
weights = np.full(t,0.5)
# dic = {
#     'Whole Blood': 0.1,
#     'Spleen': 0.1,
#     'Testis': 0.1,
#     'Esophagus - Mucosa': 0.1,
#     'Vagina': 0.1,
#     'COAD': 0.1,
#     'Lung': 0.1,
#     'Liver': 0.1,
#     'Thyroid': 0.1,
#     'BRCA': 0.1,
#     'Adipose - Subcutaneous': 0.1,
#     'Stomach': 0.1,
# }
# weights = weighting(adata,dic,weights)
weights = torch.tensor(weights,device=device)


adam = Adam({'lr': 0.002,'betas':(0.95,0.999)}) 
clipped_adam = ClippedAdam({'betas':(0.95,0.999)})
elbo = Trace_ELBO()


def model_mle(X,Y,weights):
    sigma = pyro.param('sigma',lambda:torch.tensor(np.full(n,0.5),device=device),constraint=constraints.interval(0.05,1))
    beta_y = pyro.sample('beta_y',dist.Gamma(10,1))
    beta_x = pyro.sample('beta_x',dist.Gamma(25,1))
    total = pyro.sample('total',dist.Binomial(50,weights).expand([t]).to_event(1))
    scaled_X = torch.round(X * total.unsqueeze(-1))
    with pyro.poutine.scale(scale=1), pyro.plate('data_X',t):
        c = pyro.sample('c',dist.Poisson(beta_x*sigma).expand([t,n]).to_event(1),obs=scaled_X)
    with pyro.poutine.scale(scale=1/5000), pyro.plate('data_Y',s):
        nc = pyro.sample('nc',dist.LogNormal(beta_y*sigma,0.5).expand([s,n]).to_event(1),obs=Y)


def guide_mle(X,Y,weights):
    pass

def model(X,Y,weights):
    a = torch.tensor(2.,device=device)
    b = torch.tensor(2.,device=device)
    sigma = pyro.sample('sigma',dist.Beta(a,b).expand([n]).to_event(1))
    a = torch.tensor(10.,device=device)
    b = torch.tensor(1.,device=device)
    beta_y = pyro.sample('beta_y',dist.Gamma(a,b))
    a = torch.tensor(25.,device=device)
    b = torch.tensor(1.,device=device)
    beta_x = pyro.sample('beta_x',dist.Gamma(a,b))
    a = torch.tensor(50.,device=device)
    total = pyro.sample('total',dist.Binomial(a,weights).expand([t]).to_event(1))
    scaled_X = torch.round(X * total.unsqueeze(-1))
    with pyro.poutine.scale(scale=1), pyro.plate('data_X',t):
        c = pyro.sample('c',dist.Poisson(beta_x*sigma).expand([t,n]).to_event(1),obs=scaled_X)
    with pyro.poutine.scale(scale=1/5000), pyro.plate('data_Y',s):
        nc = pyro.sample('nc',dist.LogNormal(beta_y*sigma,0.5).expand([s,n]).to_event(1),obs=Y)

def guide_map(X,Y,weights):
    sigma_param = pyro.param('sigma_param',lambda:torch.tensor(np.full(n,0.5),device=device),constraint=constraints.interval(0.05,1))
    sigma = pyro.sample('sigma',dist.Delta(sigma_param).expand([n]).to_event(1))
    a = torch.tensor(10.,device=device)
    b = torch.tensor(1.,device=device)
    beta_y = pyro.sample('beta_y',dist.Gamma(a,b))
    a = torch.tensor(25.,device=device)
    b = torch.tensor(1.,device=device)
    beta_x = pyro.sample('beta_x',dist.Gamma(a,b))
    total = pyro.sample('total',dist.Binomial(50,weights).expand([t]).to_event(1))

def guide(X,Y,weights):
    alpha = pyro.param('alpha',lambda:torch.tensor(np.full(n,2.0),device=device),constraint=constraints.positive)
    beta = pyro.param('beta',lambda:torch.tensor(np.full(n,2.0),device=device),constraint=constraints.positive)
    sigma = pyro.sample('sigma',dist.Beta(alpha,beta).expand([n]).to_event(1))
    a = torch.tensor(10.,device=device)
    b = torch.tensor(1.,device=device)
    beta_y = pyro.sample('beta_y',dist.Gamma(a,b))
    a = torch.tensor(25.,device=device)
    b = torch.tensor(1.,device=device)
    beta_x = pyro.sample('beta_x',dist.Gamma(a,b))
    total = pyro.sample('total',dist.Binomial(50,weights).expand([t]).to_event(1))
    return {'sigma':sigma,'beta_y':beta_y,'beta_x':beta_x,'total':total}

# trace = pyro.poutine.trace(model_map).get_trace(X,Y,weights)
# trace.compute_log_prob()  
# print(trace.format_shapes())
# pyro.render_model(model_map, model_args=(X,Y,weights), render_distributions=True, render_params=True, filename='model.pdf')



# train
n_steps = 5000
pyro.clear_param_store()
svi = SVI(model, guide, adam, loss=Trace_ELBO())
losses = []
for step in tqdm(range(n_steps),total=n_steps):  
    loss = svi.step(X,Y,weights)
    losses.append(loss)
plt.figure(figsize=(5, 2))
plt.plot(losses)
plt.xlabel("SVI step")
plt.ylabel("ELBO loss")
plt.savefig('elbo_loss.pdf',bbox_inches='tight')
plt.close()

with pyro.plate('samples',1000,dim=-1):
    samples = guide(X,Y,weights)
svi_sigma = samples['sigma']  # torch.Size([1000, 24290])
sigma = svi_sigma.data.cpu().numpy().mean(axis=0)
alpha = pyro.param('alpha').data.cpu().numpy()
beta = pyro.param('beta').data.cpu().numpy()
df = pd.DataFrame(index=uids,data={'mean_sigma':sigma,'alpha':alpha,'beta':beta})
with open('X.p','rb') as f:
    X = pickle.load(f)
with open('Y.p','rb') as f:
    Y = pickle.load(f)
Y_mean = Y.mean(axis=1)
X_mean = X.mean(axis=1)
df['Y_mean'] = Y_mean
df['X_mean'] = X_mean
df.to_csv('full_results.txt',sep='\t')
diagnose('full_results.txt',output_name='pyro_full_diagnosis.pdf')
sys.exit('stop')


'''evaluate'''
target = pd.read_csv('CARTargets.txt',sep='\t',index_col=0)
mapping = {e:g for g,e in target['Ensembl ID'].to_dict().items()}
target = target.loc[target['Category']=='in clinical trials',:]['Ensembl ID'].tolist()
result = pd.read_csv('mle_results.txt',sep='\t',index_col=0)
result = result.loc[target,:]
result['gene'] = result.index.map(mapping).values
result = result.sort_values(by='mean_sigma')
result.to_csv('check1.txt',sep='\t')
fig,ax = plt.subplots()
ax.bar(x=np.arange(result.shape[0]),height=result['mean_sigma'].values)
ax.set_xticks(np.arange(result.shape[0]))
ax.set_xticklabels(result['gene'].values,fontsize=1,rotation=90)
ax.set_ylabel('inferred sigma')
plt.savefig('check.pdf',bbox_inches='tight')
plt.close()


