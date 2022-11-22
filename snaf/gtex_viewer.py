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
from scipy.sparse import csr_matrix

'''
this is gtex viewer
'''

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# gtex viewer

def gtex_viewer_configuration(adata_passed_in):
    global adata
    adata = adata_passed_in

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


def gtex_visual_combine_plotly(uid,outdir='',norm=False,tumor=None):
    '''
    Generate combined plot but interactively

    :param uid: string, the uid for the splicing event
    :param norm: bool. whether normalize for the sequencing depth or not
    :param outdir: string, where the figure go into
    :param tumor: pandas dataframe, the tumor df to compare with

    Example::

        snaf.gtex_visual_combine(query='ENSG00000090339:E4.3-E4.5',norm=True,tumor=df)  
        snaf.gtex_visual_combine(query='ENSG00000112149:E7.1-E9.1',norm=True,tumor=df)

    '''
    import plotly.graph_objects as go
    query = uid
    try:
        info = adata[[query],:]
    except:
        print('{} not detected in gtex or added controls, impute as zero'.format(query))
        info = ad.AnnData(X=csr_matrix(np.full((1,adata.shape[1]),0)),obs=pd.DataFrame(data={'mean':[0]},index=[uid]),var=adata.var)  # weired , anndata 0.7.6 can not modify the X in place? anndata 0.7.2 can do that in scTriangulate
    title = query
    identifier = query.replace(':','_')
    df_gtex = pd.DataFrame(data={'value':info.X.toarray().squeeze(),'tissue':info.var['tissue'].values},index=info.var_names)
    if norm:
        df_gtex['value'] = df_gtex['value'] / (info.var['total_count'])
    if tumor is not None:
        tumor_query_value = tumor.loc[query,:].values.squeeze()
        if norm:
            tumor_total_count = tumor.sum(axis=0)
            tumor_query_value = tumor_query_value / (tumor_total_count.values.squeeze()/1e6)
        tumor_sub_df = pd.DataFrame(data={'value':tumor_query_value,'tissue':['tumor']*tumor.shape[1]},index=tumor.columns)

    expr_tumor_dict = tumor_sub_df['value'].to_dict()   # {sample:value}
    expr_tumor_dict = {sample + ',' + 'tumor': value for sample,value in expr_tumor_dict.items()}  # {sample,tumor:value}
    expr_tumor_dict = {k:v for k,v in sorted(expr_tumor_dict.items(),key=lambda x:x[1])}
    expr_gtex_df = df_gtex  # index is sample name, two columns: value and tissue
    expr_gtex_dict = {row.Index + ',' + row.tissue: row.value for row in expr_gtex_df.itertuples()}   # {sample,tissue:value}
    expr_gtex_dict = {k:v for k,v in sorted(expr_gtex_dict.items(),key=lambda x:x[1])}
    node_x = []
    node_y = []
    node_text = []
    expr_gtex_dict.update(expr_tumor_dict)
    for i,(k,v) in enumerate(expr_gtex_dict.items()):
        node_x.append(i)
        node_y.append(v)
        node_text.append(k)
    node_trace = go.Scatter(x=node_x,y=node_y,mode='markers',marker={'color':'red','size':2},text=node_text,hoverinfo='text')
    fig = go.Figure(data=[node_trace],layout=go.Layout(showlegend=False))
    fig.update_xaxes(title_text='Samples(Normal -> Tumor)')
    y_axis_title = 'Count Per Million (Normalized)' if norm else 'Raw Read Count'
    fig.update_yaxes(title_text=y_axis_title)
    fig.write_html(os.path.join(outdir,'gtex_visual_combine_plotly_norm_{}_{}.html'.format(norm,identifier)))


def gtex_visual_combine(uid,norm=False,outdir='.',figsize=(6.4,4.8),tumor=None,ylim=None):
    ''' 
    Visualize the gtex expression and tumor specificity for splicing event (combine into one plot)

    :param uid: string, the uid for the splicing event
    :param norm: bool. whether normalize for the sequencing depth or not
    :param outdir: string, where the figure go into
    :param figsize: tuple, the (width,height) of the figure
    :param tumor: pandas dataframe, the tumor df to compare with
    :param ylim: tuple, modify the ylim of the (bottom, top) of the figure

    Example::

        snaf.gtex_visual_combine(query='ENSG00000090339:E4.3-E4.5',norm=True,tumor=df)  
        snaf.gtex_visual_combine(query='ENSG00000112149:E7.1-E9.1',norm=True,tumor=df)
    '''
    query = uid
    try:
        info = adata[[query],:]
    except:
        print('{} not detected in gtex, impute as zero'.format(query))
        info = ad.AnnData(X=csr_matrix(np.full((1,adata.shape[1]),0)),obs=pd.DataFrame(data={'mean':[0]},index=[uid]),var=adata.var)  # weired , anndata 0.7.6 can not modify the X in place? anndata 0.7.2 can do that in scTriangulate
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
        sub_df.sort_values(by='value',inplace=True)
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
        ylabel = 'Normalized read counts (CPM)'
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Normal Tissues --> Tumor')
    if ylim is not None:
        ax.set_ylim(ylim)
    plt.savefig(os.path.join(outdir,'gtex_visual_combine_norm_{}_{}.pdf'.format(norm,identifier)),bbox_inches='tight')
    plt.close()
    






def gtex_visual_subplots(uid,norm=True,outdir='.'):
    ''' 
    Visualize the gtex expression and tumor specificity for splicing event (subplots)

    :param uid: string, the uid for the splicing event
    :param norm: bool. whether normalize for the sequencing depth or not
    :param outdir: string, where the figure go into

    Example::

        snaf.gtex_visual_subplots(query='ENSG00000090339:E4.3-E4.5',norm=True)
        snaf.gtex_visual_subplots(query='ENSG00000112149:E7.1-E9.1',norm=True)
    '''
    query = uid
    try:
        info = adata[[query],:]
    except:
        print('{} not detected in gtex, impute as zero'.format(query))
        info = ad.AnnData(X=csr_matrix(np.full((1,adata.shape[1]),0)),obs=pd.DataFrame(data={'mean':[0]},index=[uid]),var=adata.var)  # weired , anndata 0.7.6 can not modify the X in place? anndata 0.7.2 can do that in scTriangulate
    title = query
    identifier = query.replace(':','_')
    n_tissue = len(info.var['tissue'].unique())
    ncols = 5
    fig,axes = plt.subplots(figsize=(20,20),ncols=ncols,nrows=n_tissue//ncols+1,gridspec_kw={'wspace':0.5,'hspace':0.5})
    df_data = []
    for i,ax in enumerate(axes.flat):
        if i < n_tissue:
            tissue = info.var['tissue'].unique()[i]
            psi = info[:,info.var['tissue']==tissue].X.toarray().squeeze()
            non_zero_count = np.count_nonzero(psi)
            if norm:
                psi = psi / info[:,info.var['tissue']==tissue].var['total_count'].values
            try:
                ax.plot(np.arange(len(psi)),psi,marker='o',markersize=4,color='k',markerfacecolor='r',markeredgecolor='r',linestyle='--')
            except:
                psi = [psi]
                ax.plot(np.arange(len(psi)),psi,marker='o',markersize=4,color='k',markerfacecolor='r',markeredgecolor='r',linestyle='--')
            total_count = len(psi)
            ax.set_xticks(np.arange(len(psi)))
            ax.set_xticklabels(['s{}'.format(i) for i in np.arange(len(psi))],fontsize=4,rotation=60)
            ax.set_title('{}_count:{}/{}'.format(tissue,non_zero_count,total_count),fontsize=4)
            ax.set_ylim(bottom=-0.001,top=50)
            if norm:
                ax.set_ylabel('normalized counts')
            else:
                ax.set_ylabel('counts')
            # write the df
            df_data.append((tissue,total_count,non_zero_count,non_zero_count/total_count))
        else:
            ax.axis('off')
    df = pd.DataFrame.from_records(data=df_data,columns=['tissue','total_count','non_zero_count','ratio']).set_index(keys='tissue')
    df.to_csv(os.path.join(outdir,'tissue_specific_presence_{}.txt'.format(identifier)),sep='\t')
    fig.suptitle(title,fontsize=10)
    if norm:
        name = 'gtex_visual_subplots_norm_count_{}.pdf'.format(identifier)
    else:
        name = 'gtex_visual_subplots_count_{}.pdf'.format(identifier)
    plt.savefig(os.path.join(outdir,name),bbox_inches='tight')
    plt.close()














