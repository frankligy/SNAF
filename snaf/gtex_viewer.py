
 

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


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# gtex viewer

def gtex_viewer_configuration(adata_passed_in):
    global adata
    adata = adata_passed_in

def gtex_visual_norm_count_combined(uid,xlim=None,ylim=None,save_df=False,outdir='.'):
    '''
    This function is to visualize the normalized count (CPM) distribution,
    it should follow a half-norm distribution, this plot will help determine the choice in
    bayesian modeling for tumor specificity score

    :param uid: string, the uid for the splicing event
    :param xlim: None or tuple, whether to constrain the xlim of the histplot
    :param ylim: None or tuple, whether to constrain the ylim of the histplot
    :param save_df: boolean, whether to save a txt file reporting all the cpm values and the tissue information for debug purpose
    :param outdir: string, where the figure go into

    Example::
    
        snaf.gtex_visual_norm_count_combined(uid='ENSG00000090339:E4.3-E4.5')  
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    info = adata[[uid],:]
    scale_factor_dict = adata.var['total_count'].to_dict()
    df = pd.DataFrame(data={'value':info.X.toarray().squeeze(),'tissue':info.var['tissue'].values},index=info.var_names)
    df['value_cpm'] = df['value'].values / df.index.map(scale_factor_dict).values
    data = df['value_cpm'].values 
    fig, ax = plt.subplots()
    sns.histplot(data,bins=100,ax=ax)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.set_xlabel('Count Per Million(CPM)')
    plt.savefig(os.path.join(outdir,'hist_{}.pdf'.format(uid.replace(':','_'))),bbox_inches='tight')
    plt.close()
    if save_df:
        df.to_csv(os.path.join(outdir,'cpm_{}.txt'.format(uid.replace(':','_'))),sep='\t')
    return df

def gtex_visual_per_tissue_count(uid,total_count=10, count_cutoff=1, total=25, outdir='.'):
    '''
    This function is to visualize the random variable X, which repesents the number of samples in each
    tissue type that express this junction. It should follow a zeroinflated poisson distribution. We do not 
    consider tissues with less than total_count samples, and we use count_cutoff to represent how many count can
    be defined as "expressed" versus "non-expressed". We scaled all the number of samples to the total of "total" to make sure 
    they are comparable. Because if a tissue type has 1000 samples, another has 25 samples, this random variable is not directly
    comparable. this plot will help determine the choice in bayesian modeling for tumor specificity score

    :param uid: string, the uid for the splicing event
    :param outdir: string, where the figure go into

    Example::
    
        snaf.gtex_visual_per_tissue_count(uid='ENSG00000090339:E4.3-E4.5')  
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    x = []
    for tissue in adata.var['tissue'].unique():
        sub = adata[uid,adata.var['tissue']==tissue]
        total_count = sub.shape[1]
        if total_count >= total_count:
            c = np.count_nonzero(np.where(sub.X.toarray()<count_cutoff,0,sub.X.toarray()))
            scaled_c = round(c * (total/total_count),0)
            x.append(scaled_c)
    x = np.array(x)
    fig,ax = plt.subplots()
    try:
        sns.histplot(x,binwidth=1,ax=ax)
    except:
        sns.histplot(x,ax=ax)  # if x is with no variance at all, binwidth can not be 1
    ax.set_xlabel('Number of samples expressing this target in each tissue type')
    plt.savefig(os.path.join(outdir,'tissue_dist_{}.pdf'.format(uid.replace(':','_'))),bbox_inches='tight')
    plt.close()


def gtex_visual_combine_plotly(uid,outdir='',norm=False,tumor=None):
    '''
    Generate combined plot but interactively

    :param uid: string, the uid for the splicing event
    :param norm: bool. whether normalize for the sequencing depth or not
    :param outdir: string, where the figure go into
    :param tumor: pandas dataframe, the tumor df to compare with

    Example::

        snaf.gtex_visual_combine(uid='ENSG00000090339:E4.3-E4.5',norm=True,tumor=df)  
        snaf.gtex_visual_combine(uid='ENSG00000112149:E7.1-E9.1',norm=True,tumor=df)
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)
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


def gtex_visual_combine(uid,norm=False,outdir='.',figsize=(6.4,4.8),tumor=None,ylim=None,group_by_tissue=True):
    ''' 
    Visualize the gtex expression and tumor specificity for splicing event (combine into one plot)

    :param uid: string, the uid for the splicing event
    :param norm: bool. whether normalize for the sequencing depth or not
    :param outdir: string, where the figure go into
    :param figsize: tuple, the (width,height) of the figure
    :param tumor: pandas dataframe, the tumor df to compare with
    :param ylim: tuple, modify the ylim of the (bottom, top) of the figure
    :param group_by_tissue: bool, whether to group the normal smaples by tissue or not, default is True

    Example::

        snaf.gtex_visual_combine(uid='ENSG00000090339:E4.3-E4.5',norm=True,tumor=df)  
        snaf.gtex_visual_combine(uid='ENSG00000112149:E7.1-E9.1',norm=True,tumor=df)
    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)
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
    total_number_tissues = len(sorted_sub_df_list)
    if tumor is None:
        c_list = np.concatenate([np.array(['g']*sub_df.shape[0]) for sub_df in sorted_sub_df_list]).tolist()
    else:
        c_list_1 = np.concatenate([np.array(['g']*sub_df.shape[0]) for sub_df in sorted_sub_df_list[:-1]]).tolist()
        c_list_2 = ['r'] * sorted_sub_df_list[-1].shape[0]
        c_list = c_list_1 + c_list_2

    if not group_by_tissue:
        sorted_sub_df_list = [df,tumor_sub_df]

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
    ax.scatter(x_list,y_list,s=2,c=c_list,marker='o')
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
    plt.savefig(os.path.join(outdir,'gtex_visual_combine_norm_{}_{}_groupbytissue_{}.pdf'.format(norm,identifier,group_by_tissue)),bbox_inches='tight')
    plt.close()

    return df
    






def gtex_visual_subplots(uid,norm=True,top=100,outdir='.'):
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
            ax.set_ylim(bottom=-0.001,top=top)
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









