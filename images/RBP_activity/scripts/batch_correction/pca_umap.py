#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os
import sys
import math
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import umap
import matplotlib.patches as mpatches
from colors import *



## PCA and UMAP analysis
def scikit_learn_pca(X,p):
    from sklearn.decomposition import PCA
    n_components = p
    model = PCA(n_components=n_components)
    scoring = model.fit_transform(X)   # [n_samples,n_components]
    loading = model.components_.T    # [n_events,n_components]
    variance_ratio = model.explained_variance_ratio_  # [n_components,]
    return scoring,loading,variance_ratio


def draw_pca_umap(adata,name,cat=['batch','is_control'],outdir='./'):
    data = adata.X.T
    scoring,loading,variance_ratio = scikit_learn_pca(data,p=100)
    reducer = umap.UMAP(random_state=42,min_dist=0.1)
    embedding = reducer.fit_transform(scoring)

    fig,axes = plt.subplots(nrows=2,ncols=2,gridspec_kw={'wspace':0.5,'hspace':0.5},figsize=(20,20))

    ## add col color bar
    
    for column in cat:
        color_dict = colors_for_set(adata.var[column].unique().tolist())
        adata.var['{}_color'.format(column)] = adata.var[column].map(color_dict).values

    # pca_colorby_batch
    ax = axes[0,0]
    ax.set_title('pca_colorby_batch')
    ax.scatter(scoring[:,0],scoring[:,1],c=adata.var['batch_color'].values,s=6)
    ax.set_xlabel('PC1 ({})'.format(variance_ratio[0]))
    ax.set_ylabel('PC2 ({})'.format(variance_ratio[1]))
    ax.legend(handles=[mpatches.Patch(color=i) for i in adata.var['batch_color'].unique()],labels=[i for i in adata.var['batch'].unique()],
                loc='upper left',bbox_to_anchor=(1,1))

    # pca_colorby_iscontrol
    ax = axes[0,1]
    ax.set_title('pca_colorby_iscontrol')
    ax.scatter(scoring[:,0],scoring[:,1],c=adata.var['is_control_color'].values,s=6)
    ax.set_xlabel('PC1 ({})'.format(variance_ratio[0]))
    ax.set_ylabel('PC2 ({})'.format(variance_ratio[1]))
    ax.legend(handles=[mpatches.Patch(color=i) for i in adata.var['is_control_color'].unique()],labels=[i for i in adata.var['is_control'].unique()],
                loc='upper left',bbox_to_anchor=(1,1))


    # umap_colorby_batch
    ax = axes[1,0]
    ax.set_title('umap_colorby_batch')
    ax.scatter(embedding[:,0],embedding[:,1],c=adata.var['batch_color'].values,s=6)
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')


    # umap_colorby_iscontrol
    ax = axes[1,1]
    ax.set_title('umap_colorby_iscontrol')
    ax.scatter(embedding[:,0],embedding[:,1],c=adata.var['is_control_color'].values,s=6)
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    plt.savefig(os.path.join(outdir,'pca_umap_{}.pdf'.format(name)),bbox_inches='tight')
    plt.close()


def draw_all_pca_umap(adata,name,column,outdir='./',plot_umap=True,ncols=5,layer=None):
    if layer is None:
        data = adata.X.T
    else:
        data = adata.layers[layer].T
    scoring,loading,variance_ratio = scikit_learn_pca(data,p=100)  

    batches = adata.var[column].unique()
    n_batch = len(adata.var[column].unique())
    n_plot = n_batch + 1 # add the whole plot
    nrows = math.ceil(n_plot / ncols)
    fig,axes = plt.subplots(nrows=nrows,ncols=ncols,gridspec_kw={'wspace':0.5,'hspace':0.5},figsize=(20,20))
    flattened_axes = axes.flatten()

    color_dict = colors_for_set(adata.var[column].unique().tolist())
    adata.var['{}_color'.format(column)] = adata.var[column].map(color_dict).values 

    for i,ax in enumerate(flattened_axes):
        if i == 0:
            ax.set_title('pca_colorby_batch')
            ax.scatter(scoring[:,0],scoring[:,1],c=adata.var['{}_color'.format(column)].values,s=3)
            ax.set_xlabel('PC1 ({})'.format(variance_ratio[0]))
            ax.set_ylabel('PC2 ({})'.format(variance_ratio[1]))
        elif i <= n_batch: 
            batch = batches[i-1]  # plot1, means batch 0, since plot 0 is the overview
            ax.set_title('{}'.format(batch))  
            ax.scatter(scoring[:,0],scoring[:,1],c= ['r' if item == batch else 'grey' for item in adata.var[column]],s=3)
        else:
            ax.axis('off')
    
    plt.savefig(os.path.join(outdir,'{}_pca_{}.pdf'.format(name,column)),bbox_inches='tight')
    plt.close()

    if plot_umap:
        reducer = umap.UMAP(random_state=42,min_dist=0.1)
        embedding = reducer.fit_transform(scoring) 
        fig,axes = plt.subplots(nrows=nrows,ncols=ncols,gridspec_kw={'wspace':0.5,'hspace':0.5},figsize=(20,20))
        flattened_axes = axes.flatten() 
        for i,ax in enumerate(flattened_axes):
            if i == 0:
                ax.set_title('umap_colorby_batch')
                ax.scatter(embedding[:,0],embedding[:,1],c=adata.var['{}_color'.format(column)].values,s=3)
            elif i<= n_batch: 
                batch = batches[i-1]  # plot1, means batch 0, since plot 0 is the overview
                ax.set_title('{}'.format(batch))  
                ax.scatter(embedding[:,0],embedding[:,1],c= ['r' if item == batch else 'grey' for item in adata.var[column]],s=3)
            else:
                ax.axis('off')
        
        plt.savefig(os.path.join(outdir,'{}_umap_{}.pdf'.format(name,column)),bbox_inches='tight')
        plt.close()






            






