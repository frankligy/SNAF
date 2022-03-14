#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt
import anndata as ad
import seaborn as sns

'''
this is gtex viewer
'''


# gtex viewer

def gtex_viewer_configuration(gtex_db):
    global adata
    adata = ad.read_h5ad(gtex_db)

def gtex_visual_norm_count_combined(query,out_folder='.'):
    data = adata[query,:].X.toarray().squeeze() / adata.var['total_count'].values
    sns.histplot(data,binwidth=0.01)
    plt.savefig(os.path.join(out_folder,'hist_{}.pdf'.format(query.replace(':','_'))),bbox_inches='tight')
    plt.close()

def gtex_visual_per_tissue_count(query,out_folder='.'):
    per_tissue_count = []
    for tissue in adata.var['tissue'].unique():
        sub = adata[query,adata.var['tissue']==tissue]
        c = np.count_nonzero(sub.X.toarray())
        per_tissue_count.append(c)
    sns.histplot(np.array(per_tissue_count),binwidth=1)
    plt.savefig(os.path.join(out_folder,'poisson_{}.pdf'.format(query.replace(':','_'))),bbox_inches='tight')
    plt.close()


def gtex_visual_combine(query,norm=False,outdir='.',figsize=(6.4,4.8),tumor=None):
    ''' 
    Example:
    snaf.gtex_visual_combine(query='ENSG00000090339:E4.3-E4.5',norm=True,tumor=df)  
    snaf.gtex_visual_combine(query='ENSG00000112149:E7.1-E9.1',norm=True,tumor=df)
    '''
    try:
        info = adata[[query],:]
    except:
        print('{} not detected in gtex, impute as zero'.format(query))
        info = adata[['ENSG00000090339:E4.3-E4.5'],:]
        info.X = np.full((1,info.shape[1]),0)
    title = query
    identifier = query.replace(':','_')
    df = pd.DataFrame(data={'value':info.X.toarray().squeeze(),'tissue':info.var['tissue'].values},index=info.var_names)
    if norm:
        df['value'] = df['value'] / (info.var['total_count'])
    tmp = []
    for tissue,sub_df in df.groupby(by='tissue'):
        tmp.append((sub_df,sub_df['value'].mean()))
    sorted_sub_df_list = list(list(zip(*sorted(tmp,key=lambda x:x[1])))[0])
    if tumor is not None:
        tumor_query_value = tumor.loc[query,:].values.squeeze()
        if norm:
            tumor_total_count = tumor.sum(axis=0)
            tumor_query_value = tumor_query_value / (tumor_total_count.values.squeeze()/1e6)
        tumor_sub_df = pd.DataFrame(data={'value':tumor_query_value,'tissue':['tumor']*tumor.shape[1]},index=tumor.columns)
        sorted_sub_df_list.append(tumor_sub_df)
    fig,ax = plt.subplots(figsize=figsize)
    x = 0
    x_list = []
    y_list = []
    v_delimiter = [0]
    xticklabel = []
    for i,sub_df in enumerate(sorted_sub_df_list):
        n = sub_df.shape[0]
        xticklabel.append(sub_df['tissue'].iloc[0])
        for j,v in enumerate(sub_df['value']):
            x_list.append(x)
            y_list.append(v)
            x += 1
            if j == n-1:
                v_delimiter.append(x)
                x += 1
    ax.plot(x_list,y_list,marker='o',linestyle='',markerfacecolor='r',markeredgewidth=0.1,color='k',markersize=2)
    for v in v_delimiter[1:-1]:
        ax.axvline(v,linestyle='--',linewidth=0.5)
    xtick = [(v + v_delimiter[i+1])/2 for i,v in enumerate(v_delimiter[:-1])]
    ax.set_xticks(xtick)
    ax.set_xticklabels(xticklabel,rotation=90,fontsize=1)
    ax.set_title(title)
    ylabel = 'Raw read counts'
    if norm:
        ylabel = 'Normalized read counts (TPM)'
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Normal Tissues --> Tumor')
    plt.savefig(os.path.join(outdir,'gtex_visual_combine_norm_{}_{}.pdf'.format(norm,identifier)),bbox_inches='tight')
    plt.close()
    






def gtex_visual_subplots(query,norm=True,outdir='.'):
    ''' 
    Example:
    snaf.gtex_visual_subplots(query='ENSG00000090339:E4.3-E4.5',norm=True)
    snaf.gtex_visual_subplots(query='ENSG00000112149:E7.1-E9.1',norm=True)
    '''
    try:
        info = adata[[query],:]
    except:
        print('{} not detected in gtex, impute as zero'.format(query))
        info = adata[['ENSG00000090339:E4.3-E4.5'],:]
        info.X = np.full((1,info.shape[1]),0)
    title = query
    identifier = query.replace(':','_')
    n_tissue = len(info.var['tissue'].unique())
    ncols = 5
    fig,axes = plt.subplots(figsize=(20,20),ncols=ncols,nrows=n_tissue//ncols+1,gridspec_kw={'wspace':0.5,'hspace':0.5})
    for i,ax in enumerate(axes.flat):
        if i < n_tissue:
            tissue = info.var['tissue'].unique()[i]
            psi = info[:,info.var['tissue']==tissue].X.toarray().squeeze()
            if norm:
                psi = psi / info[:,info.var['tissue']==tissue].var['total_count'].values
            ax.plot(np.arange(len(psi)),psi,marker='o',markersize=4,color='k',markerfacecolor='r',markeredgecolor='r',linestyle='--')
            ax.set_xticks(np.arange(len(psi)))
            ax.set_xticklabels(['s{}'.format(i) for i in np.arange(len(psi))],fontsize=4,rotation=60)
            ax.set_title(tissue,fontsize=8)
            ax.set_ylim(bottom=-0.001,top=50)
            if norm:
                ax.set_ylabel('normalized counts')
            else:
                ax.set_ylabel('counts')
        else:
            ax.axis('off')
    fig.suptitle(title,fontsize=10)
    if norm:
        name = 'gtex_visual_norm_count_{}.pdf'.format(identifier)
    else:
        name = 'gtex_visual_count_{}.pdf'.format(identifier)
    plt.savefig(os.path.join(outdir,name),bbox_inches='tight')
    plt.close()














