#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import statsmodels.stats.multitest as ssm
from scipy.stats import mannwhitneyu,pearsonr,spearmanr
from tqdm import tqdm


'''
this script contains survival analysis, mutation analysis
'''

def mutation_analysis(mode,burden,mutation,output,n_sample_cutoff=10,gene_column='gene',genes_to_plot=None):
    # sample is in the index of mutation
    # burden is a series, make sure the index are consistent
    if mode == 'compute':
        burden = burden.loc[burden.index.isin(mutation.index)]
        with open('snaf_mutation_tmp','w') as f:
            f.write('mutation_gene\tn_samples\tpval\n')
            for gene in tqdm(mutation['gene'].unique()):
                yes_samples = mutation.loc[mutation[gene_column] == gene,:].index.unique().tolist()
                if len(yes_samples) < n_sample_cutoff:
                    continue
                burden_df = burden.to_frame()
                burden_df.columns = ['burden']
                burden_df['mutation_{}'.format(gene)] = [True if sample in set(yes_samples) else False for sample in burden_df.index]
                x = burden_df.loc[burden_df['mutation_{}'.format(gene)],'burden'].values
                y = burden_df.loc[~(burden_df['mutation_{}'.format(gene)]),'burden'].values
                u,p = mannwhitneyu(x=x,y=y)
                f.write('{}\t{}\t{}\n'.format(gene,len(yes_samples),p))
            asso = pd.read_csv('snaf_mutation_tmp',sep='\t',index_col=0)
            results = ssm.multipletests(asso['pval'].values,alpha=0.05,method='fdr_bh')
            asso['adjp'] = results[1]
            asso.to_csv(output,sep='\t')
            os.remove('snaf_mutation_tmp')
    elif mode == 'plot':
        for gene in genes_to_plot:
            yes_samples = mutation.loc[mutation[gene_column] == gene,:].index.unique().tolist()
            burden_df = burden.to_frame()
            burden_df.columns = ['burden']
            burden_df['mutation_{}'.format(gene)] = [True if sample in set(yes_samples) else False for sample in burden_df.index]
            x = burden_df.loc[burden_df['mutation_{}'.format(gene)],'burden'].values
            y = burden_df.loc[~(burden_df['mutation_{}'.format(gene)]),'burden'].values
            u,p = mannwhitneyu(x=x,y=y)
            fig,ax = plt.subplots()
            sns.boxplot(data=burden_df,x='mutation_{}'.format(gene),y='burden',ax=ax,width=0.5)
            ax.text(x=0.5,y=0.5,s='mannwhitney p={}'.format(round(p,4)),weight='bold')
            plt.savefig('{}_{}.pdf'.format(output,gene),bbox_inches='tight')
            plt.close()



def survival_analysis(burden,survival,n,stratification_plot,survival_plot,
                      survival_duration='OS.time',survival_event='OS'): 
    # burden is a seires, survival is a df, make sure the index are the consistent
    # sample is in the index of survival
    burden = burden.loc[burden.index.isin(survival.index)]
    if n == 4:
        quantiles = burden.quantile([0.25,0.5,0.75]).values
        iqr = quantiles[2] - quantiles[0]
        upper_bound = quantiles[2] + 1.5*iqr
        lower_bound = quantiles[0] - 1.5*iqr
        identity_col = []
        for item in burden:
            if item > upper_bound:
                identity_col.append('outlier')
            elif item > quantiles[2] and item <= upper_bound:
                identity_col.append('high')
            elif item > quantiles[0] and item <= quantiles[2]:
                identity_col.append('medium')
            elif item >= lower_bound and item <= quantiles[0]:
                identity_col.append('low')
            elif item < lower_bound:
                identity_col.append('outlier')
    elif n == 3:
        quantiles = burden.quantile([0.33,0.66]).values
        identity_col = []
        for item in burden:
            if item > quantiles[1]:
                identity_col.append('high')
            elif item > quantiles[0] and item <= quantiles[1]:
                identity_col.append('medium')
            elif item <= quantiles[0]:
                identity_col.append('low')
    elif n == 2:
        quantiles = burden.quantile([0.5]).values[0]  # a scalar
        identity_col = []
        for item in burden:
            if item > quantiles:
                identity_col.append('high')
            else:
                identity_col.append('low')
    burden_encode = pd.Series(index=burden.index,data=identity_col)
    be_vc = burden_encode.value_counts()
    # plot stratification
    sns.boxplot(data=burden.values);plt.savefig(stratification_plot,bbox_inches='tight');plt.close()
    high_group = burden_encode.loc[burden_encode=='high'].index.tolist()
    low_group = burden_encode.loc[burden_encode=='low'].index.tolist()
    high_os = survival.loc[high_group,[survival_duration,survival_event]]
    low_os = survival.loc[low_group,[survival_duration,survival_event]]
    # plot KM plot and logrank test
    fig,ax = plt.subplots()
    ax.set_ylim(0,1)
    for df in [low_os,high_os]:
        kmf = KaplanMeierFitter()
        kmf.fit(df[survival_duration],df[survival_event])
        kmf.plot_survival_function(ax=ax,ci_show=False,at_risk_counts=False)
    current_handles,current_labels = ax.get_legend_handles_labels()
    new_labels = ['low_burden','high_burden']
    ax.legend(current_handles,new_labels)
    results = logrank_test(low_os[survival_duration],high_os[survival_duration],low_os[survival_event],high_os[survival_event])
    ax.text(x=1000,y=0.05,s='Log-rank test: p-value is {:.2f}'.format(results.p_value),weight='bold')
    plt.savefig(survival_plot,bbox_inches='tight');plt.close()
    return burden,burden_encode,be_vc

def stage0_compatible_results(jcmq,outdir='.',name_burden='burden_stage0.txt',name_frequency='frequency_stage0.txt'):
    # for burden
    burden_stage0 = jcmq.cond_df.astype('int8')
    burden_stage0.loc['burden'] = burden_stage0.sum(axis=0).values
    burden_stage0['mean'] = burden_stage0.mean(axis=1).values
    burden_stage0.to_csv(os.path.join(outdir,name_burden),sep='\t')
    # for frequency
    tuple_list = []
    for index,series in jcmq.cond_df.iterrows():
        samples = series.loc[series].index.tolist()
        tuple_ = (index,samples,len(samples))
        tuple_list.append(tuple_)
    df = pd.DataFrame.from_records(tuple_list,columns=['junction','samples','n_sample']).set_index(keys='junction').sort_values(by='n_sample',ascending=False)
    df.to_csv(os.path.join(outdir,name_frequency),sep='\t')





def reformat_frequency_table(df):
    # the index has to be comma separated string
    from ast import literal_eval
    df['samples'] = [literal_eval(item) for item in df['samples']]
    sequence_io = []
    for row in df.itertuples():
        for item in row.samples:
            sequence_io.append((row.Index,item,1))
    result = pd.DataFrame.from_records(sequence_io,columns=['id','sample','value'])
    result = result.groupby(by=['id','sample'])['value'].mean().unstack(fill_value=0)
    return result
        

    