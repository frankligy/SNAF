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
from copy import deepcopy
import multiprocessing as mp

'''
this script contains survival analysis, mutation analysis
'''

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

def plot_umap_neoantigen(df_path,outdir):
    df = pd.read_csv(df_path,sep='\t',index_col=0)
    fig,ax = plt.subplots()
    ax.scatter(df['umap_x'],df['umap_y'],color=['red' if item == 'high' else 'blue' for item in df['identity']],s=3)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([],[],marker='o',linestyle='',color=i) for i in ['red','blue']],labels=['Shared Neoantigen','Unique Neoantigen'],bbox_to_anchor=(1,1),loc='upper left',frameon=False)
    ax.set_xlabel('umap_x')
    ax.set_ylabel('umap_y')
    plt.savefig(os.path.join(outdir,'mer_umap.pdf'),bbox_inches='tight')
    plt.close()

def survival_regression(freq,remove_quote,rename_func,survival,pea,outdir='.',cores=None,mode='binary',survival_duration='OS.time',survival_event='OS',n_effective_obs=3):
    '''
    conduct cox regression to identify neoantigens whose parental junctions expression is significantly associated with survival

    :param freq: string, path to the T antigen generated frequency dataframe
    :param remove_quote: bool, whether to remove the quote of samples column for freq df or not
    :param rename_func: function, a function to rename columns in freq and event annotation file to be consistent with survival file
    :param survival: string, the path to the survival dataframe
    :param pea: string, the path to the altanalyze generated event annotation file
    :param outdir: string, the output directory
    :param cores: None or int, how many cores to use for computation
    :param mode: binary or psi, if binary, we test whether having a neoantigen can predict survival, psi means whether its parental junction's psi value is associated with survival
    :param survival_duration: string, the column name for duration in survival dataframe
    :param survival_event: string, the column name for event in survival dataframe
    :param n_effective_obs: int, default is 3, at least three obs whose event=1, otherwise all obs are right-sensored

    Examples::

        snaf.downstream.survival_regression(freq='result/frequency_stage3_verbosity1_uid_gene_symbol_coord_mean_mle.txt',remove_quote=True,
                                            rename_func=lambda x:'-'.join(x.split('-')[:4]),survival='TCGA-SKCM.survival.tsv',
                                            pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',outdir='result/survival')
    '''
    # reformat psi matrix for further extraction
    ea = pd.read_csv(pea,sep='\t',index_col='UID').iloc[:,10:]
    ea.index = [':'.join(item.split('|')[0].split(':')[1:]) for item in ea.index]
    ea = ea.loc[np.logical_not(ea.index.duplicated()).tolist(),:]
    ea.rename(columns=rename_func,inplace=True)
    # format freq matrix
    freq = pd.read_csv(freq,index_col=0,sep='\t')
    from ast import literal_eval
    if remove_quote:
        freq['samples'] = [literal_eval(item) for item in freq['samples']]
    # read in survival data
    survival = pd.read_csv(survival,index_col=0,sep='\t')
    survival = survival.loc[:,[survival_duration,survival_event]]
    # start to calculate association for each neoantigen 
    # before spawn, let's rename all the samples
    freq['samples'] = [[rename_func(s) for s in ss] for ss in freq['samples']]
    if cores is None:
        cores = mp.cpu_count()
    pool = mp.Pool(processes=cores)
    print('{} subprocesses have been spawned'.format(cores))
    sub_dfs = split_df_to_chunks(freq,cores=cores)
    if mode == 'binary':
        r = [pool.apply_async(func=survival_regression_binary_atomic,args=(sub_df,ea,survival,survival_duration,survival_event,n_effective_obs,)) for sub_df in sub_dfs]  
        pool.close()
        pool.join()
        df_data = []
        for collect in r:
            result = collect.get()
            df_data.extend(result)
    elif mode == 'psi':
        r = [pool.apply_async(func=survival_regression_psi_atomic,args=(sub_df,ea,survival,survival_duration,survival_event,n_effective_obs,)) for sub_df in sub_dfs]  
        pool.close()
        pool.join()
        df_data = []
        for collect in r:
            result = collect.get()
            df_data.extend(result)    
    final_df = pd.DataFrame.from_records(data=df_data,columns=['uid','pep','junc','n_valid_sample','wald_z_score','wald_p_value'])
    final_df.set_index(keys='uid',inplace=True)
    final_df.to_csv(os.path.join(outdir,'survival_regression_{}_final_results.txt'.format(mode)),sep='\t')

def split_df_to_chunks(df,cores=None):
    df_index = np.arange(df.shape[0])
    if cores is None:
        cores = mp.cpu_count()
    sub_indices = np.array_split(df_index,cores)
    sub_dfs = [df.iloc[sub_index,:] for sub_index in sub_indices]
    return sub_dfs

def survival_regression_binary_atomic(freq,ea,survival,survival_duration,survival_event,n_effective_obs):
    df_data = []
    from lifelines import CoxPHFitter
    for uid,ss in tqdm(zip(freq.index,freq['samples']),total=freq.shape[0]):
        pep,junc = uid.split(',')
        tmp = ea.loc[junc,:].to_frame()
        tmp.columns = ['psi']
        # binary, as psi is variable, we just want to know if they are present or not
        tmp.loc[ss,['psi']] = 1
        not_ss = list(set(tmp.index).difference(set(ss)))
        tmp.loc[not_ss,['psi']] = 0
        operation_df = tmp.join(other=survival,how='inner')
        operation_df = operation_df.dropna()  # actually unneccessary, because if a neojunction can give rise to neoantigen, so neojunction psi should be above an appreciable level
        n_valid_sample = len(ss)
        try:
            assert operation_df.loc[operation_df[survival_event]==1,:].shape[0] >= n_effective_obs
        except AssertionError:
            wald_z_score, wald_p_value = 'n_effective_obs<3','n_effective_obs<3'
            df_data.append((uid,pep,junc,n_valid_sample,wald_z_score,wald_p_value))
        else:
            cph = CoxPHFitter()
            try:
                cph.fit(operation_df,duration_col=survival_duration,event_col=survival_event)
                summary = cph.summary
                wald_z_score, wald_p_value = summary.loc['psi',['z','p']].tolist()
            except:
                wald_z_score, wald_p_value = 'convergence_error','convergence_error'
            df_data.append((uid,pep,junc,n_valid_sample,wald_z_score,wald_p_value))
    return df_data    

def survival_regression_psi_atomic(freq,ea,survival,survival_duration,survival_event,n_effective_obs):
    df_data = []
    from lifelines import CoxPHFitter
    for uid,ss in tqdm(zip(freq.index,freq['samples']),total=freq.shape[0]):
        pep,junc = uid.split(',')
        tmp = ea.loc[junc,ss].to_frame()
        tmp.columns = ['psi']
        operation_df = tmp.join(other=survival,how='inner')
        if uid == 'RETDFKMKF,ENSG00000167291:E38.6-E39.1':
            operation_df.to_csv('tmp_test_old.txt',sep='\t');sys.exit('stop')
        operation_df = operation_df.dropna()  # actually unneccessary, because if a neojunction can give rise to neoantigen, so neojunction psi should be above an appreciable level
        n_valid_sample = operation_df.shape[0]
        try:
            assert operation_df.loc[operation_df[survival_event]==1,:].shape[0] >= n_effective_obs
        except AssertionError:
            wald_z_score, wald_p_value = 'n_effective_obs<3','n_effective_obs<3'
            df_data.append((uid,pep,junc,n_valid_sample,wald_z_score,wald_p_value))
        else:
            cph = CoxPHFitter()
            try:
                cph.fit(operation_df,duration_col=survival_duration,event_col=survival_event)
                summary = cph.summary
                wald_z_score, wald_p_value = summary.loc['psi',['z','p']].tolist()
            except:
                wald_z_score, wald_p_value = 'convergence_error','convergence_error'
            df_data.append((uid,pep,junc,n_valid_sample,wald_z_score,wald_p_value))
    return df_data




def prepare_DEG_analysis(burden_path,patient_strat_path,exp_path,outdir='.',encoding={'low':'1','high':'2'}):
    '''
    This function generate groups and comps file for performing AltAnalyze Differential Gene Analysis using empirical Bayes

    :param burden_path: string, the path to the burden file generated by T pipeline
    :param patient_strat_path: string, the path to the patient stratification file returned by survival_analysis function
    :param exp_path: string, the gene expression file (exp.original.txt) generated by step1 AltAnalyze results, should be in ExpressionInput
    :param outdir: string, the folder for generating comps and groups files
    :param encoding: dict, the encoding to specificy how to transform identity string and what identities will be used for the DEG analysis


    Examples::

        snaf.downstream.prepare_DEG_analysis('result/burden_stage3.txt','result/survival/burden3_patient_high_low_group.txt',
                                             '/data/salomonis2/NCI-R01/TCGA-SKCM/TCGA_SKCM_hg38/ExpressionInput/exp.TCGA-SKCM.txt',
                                             outdir='result/survival',encoding={'low':'1','high':'2'})


    '''
    burden3 = pd.read_csv(burden_path,sep='\t',index_col=0)
    old_to_new = {}
    new_to_old = {}
    for x in burden3.columns:
        new_x = '-'.join(x.split('-')[:4])
        old_to_new[x] = new_x
        new_to_old[new_x] = x
    strat = pd.read_csv(patient_strat_path,sep='\t',index_col=0)
    ident_dict = encoding
    with open(os.path.join(outdir,'comps.txt'),'w') as f:
        f.write('{}\t{}\n'.format(list(ident_dict.values())[0],list(ident_dict.values())[1]))
    with open(os.path.join(outdir,'groups.txt'),'w') as f:
        for row in strat.itertuples():
            if row.identity == list(ident_dict.keys())[0] or row.identity == list(ident_dict.keys())[0]:
                string = new_to_old[row.Index] + '\t' + ident_dict[row.identity] + '\t' + row.identity
                f.write('{}\n'.format(string))

    
def visualize_DEG_result(result_path,mode,outdir='.',genes_to_highlight=[],up_cutoff=0.8,down_cutoff=-0.8,xlims=(-2,2),ylims=None):
    '''
    Generate volcano plot for the AltAnalyze DEG results

    :param result_path: string, the path to the AltAnalyze DEG results
    :param mode: string, either interactive or static
    :param outdir: string, the folder where the plots will go into
    :param genes_to_highlight: list, the gene names to highlight in the plot, only for static mode
    :param up_cutoff: float, the up-regulated cutoff for logfoldchange, only for static mode
    :param down_cutoff: float, the down-regulated cutoff for logfoldchagne, only for static mode
    :param xlims: None or tuple, force the xlims of the volcano plot, only for static mode
    :param ylims: None or tuple, force the ylims of the volcano plot, only for static mode

    Examples:

        snaf.downstream.visualize_DEG_result('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result/survival/DEGs-LogFold_0.58_adjp/GE.low_vs_high.txt',
                                             mode='static',outdir='.',genes_to_highlight=['LST1','HCST','IL32','CD3D','S100A8','MZB1','IGLC4','ADAM10','ARFGEF2','MIB1','KIF3B','TNPO1','PTPN11','ANKRD52','TGFBR1'])

    '''
    
    df = pd.read_csv(result_path,sep='\t',index_col=0).loc[:,['LogFold','adjp','Symbol']]
    df['minus_log10_adjp'] = np.negative(np.log10(df['adjp'].values))
    df.drop_duplicates(subset=['Symbol'],inplace=True)
    problematics = list(set(genes_to_highlight).difference(set(df['Symbol'].values)))
    print('{} are not in DEG results'.format(problematics))
    if mode == 'interactive':
        import plotly.graph_objects as go
        node_x = []
        node_y = []
        node_text = []
        for row in df.itertuples():
            node_x.append(float(row.LogFold))
            node_y.append(float(row.minus_log10_adjp))
            node_text.append(row.Symbol)
        node_trace = go.Scatter(name='nodes',x=node_x,y=node_y,mode='markers',hoverinfo='text',text=node_text,marker={'color':'blue','size':5})
        fig_layout = go.Layout(showlegend=False,title='DEG',xaxis=dict(title_text='logfold'),yaxis=dict(title_text='-log10_adjp'))
        fig = go.Figure(data=[node_trace],layout=fig_layout)
        fig.write_html(os.path.join(outdir,'volcano_plot_DEG_interactive.html'),include_plotlyjs='cdn')
    elif mode == 'static':
        from adjustText import adjust_text
        from matplotlib.patches import ArrowStyle
        node_x = []
        node_y = []
        node_text = []
        node_color = []
        for row in df.itertuples():
            node_x.append(float(row.LogFold))
            node_y.append(float(row.minus_log10_adjp))
            node_text.append(row.Symbol)
            if float(row.LogFold) > up_cutoff:
                node_color.append('#F22727')
            elif float(row.LogFold) < down_cutoff:
                node_color.append('#41D2F2')
            else:
                node_color.append('#CFCFD0')
        fig,ax = plt.subplots()
        ax.scatter(node_x,node_y,s=2,c=node_color)  
        ax.spines['left'].set_position('center')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel('logfold')
        ax.set_ylabel('-log10(adjp)')
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)
        geneset = set(genes_to_highlight)
        texts = [ax.text(node_x[i],node_y[i],node_text[i],ha='center',va='center',fontsize=5) for i in range(len(node_text)) if node_text[i] in geneset]
        adjust_text(texts,arrowprops=dict(arrowstyle='-',color='#CFCFD0'))
        plt.savefig(os.path.join(outdir,'volcano_plot_DEG_static.pdf'),bbox_inches='tight');plt.close()


def prepare_GO_analysis(result_path,lc_cutoff=0.5,adjp_cutoff=0.05,n=None,sortby='adjp',ascending=True,outdir='.'):
    '''
    Generate gene list for GO analysis

    :param result_path: string, the path to the AltAnalyze generated DEG result
    :param lc_cutoff, float, the cutoff above which will be considered desirable genes
    :param adjp_cutoff, float, the cutoff below which will be considered desirable genes
    :param n: None or int, the number of genes to contrain for GO analysis
    :param sortby: string, 'adjp' or 'LogFold'
    :param ascending: bool, sort in ascending order or not
    :param outdir: string, the directory to write outputs

    Examples::

        snaf.downstream.prepare_GO_analysis('result/survival/DEGs-LogFold_0.58_adjp/GE.low_vs_high.txt',outdir='result/survival')


    '''
    df = pd.read_csv(result_path,sep='\t',index_col=0).loc[:,['LogFold','adjp','Symbol']]
    df.drop_duplicates(subset=['Symbol'],inplace=True)
    valid_df = df.loc[(df['LogFold']>lc_cutoff) & (df['adjp']<adjp_cutoff),:]
    print('Total genes are {}, set n as {}, sort by {}, ascending {}'.format(valid_df.shape[0],n,sortby, ascending))
    valid_df = valid_df.sort_values(by=sortby,axis=0,ascending=ascending)
    valid_df.to_csv(os.path.join(outdir,'gene_dataframe.txt'),sep='\t')
    gene_list = valid_df['Symbol'].tolist()
    if n is None:
        n = len(gene_list)
    else:
        if n > len(gene_list):
            n = len(gene_list)
    with open(os.path.join(outdir,'gene_list.txt'),'w') as f:
        c = 0
        for i,gene in enumerate(gene_list):
            if c < n:
                f.write('{}\n'.format(gene))
                c += 1

def visualize_GO_result(path_list,category_list=['Gene-Set Name','Ontology Name(Ontology-ID)'],mode='interactive',outdir='',ontology_to_highlight={},
                        xlims=None,ylims=None):
    '''
    Visualize the GO analysis rssults

    :param path_list: list, each element points to the GO-Elite result txt file
    :param category_list, list, each element points to the column name in each result file that represent GO terms
    :param mode: string, interactive or static
    :param outdir: string, the output directory
    :param ontology_to_highlight: dict, key is the GO terms in the result, values are the corresponding text you want them to be displayed on the plot
    :param xlims: None or tuple.
    :param ylims: None or tuple.

    Examples::

        snaf.downstream.visualize_GO_result(path_list=['result/survival/GO_Elite_result_BioMarkers/GO-Elite_results/pruned-results_z-score_elite.txt','result/survival/GO_Elite_result_GeneOntology/GO-Elite_results/pruned-results_z-score_elite.txt'],
                                    mode='static',ontology_to_highlight={'Placenta T cells2 (PMID32214235 top 100)':'Placenta T cell','Fetal Liver DC2 (PMC6861135 markers)':'Fetal Liver DC','complement activation(GO:0006956)':'Complement Activation'},ylims=(-0.05,0.1))

    '''
    dfs = []
    for path,column in zip(path_list,category_list):
        df = pd.read_csv(path,sep='\t',index_col=0).loc[:,[column,'Z Score','AdjustedP']]
        df.rename(columns={column:'gene_set_name','Z Score':'z_score'},inplace=True)
        dfs.append(df)
    df_all = pd.concat(dfs,axis=0)
    if mode == 'interactive':
        import plotly.graph_objects as go
        node_x = []
        node_y = []
        node_text = []
        for row in df_all.itertuples():
            node_x.append(float(row.z_score))
            node_y.append(float(row.AdjustedP))
            node_text.append(row.gene_set_name)
        node_trace = go.Scatter(name='nodes',x=node_x,y=node_y,mode='markers',hoverinfo='text',text=node_text,marker={'color':'blue','size':5})
        fig_layout = go.Layout(showlegend=False,title='DEG',xaxis=dict(title_text='z_score',range=xlims),yaxis=dict(title_text='FDR_p',range=ylims))
        fig = go.Figure(data=[node_trace],layout=fig_layout)
        fig.write_html(os.path.join(outdir,'GO_plot_interactive.html'),include_plotlyjs='cdn')
    elif mode == 'static':
        from adjustText import adjust_text
        from matplotlib.patches import ArrowStyle
        node_x = []
        node_y = []
        node_text = []
        for row in df_all.itertuples():
            node_x.append(float(row.z_score))
            node_y.append(float(row.AdjustedP))
            node_text.append(row.gene_set_name)
        fig,ax = plt.subplots()
        ax.scatter(node_x,node_y,s=2,c='#3D52A3')  
        ax.spines['left'].set_position(('axes',0.3))
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel('GO-Elite Z-Score')
        ax.set_ylabel('FDR_p')
        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)
        pathway_set = set(ontology_to_highlight.keys())
        texts = [ax.text(node_x[i],node_y[i],ontology_to_highlight[node_text[i]],ha='center',va='center',fontsize=3) for i in range(len(node_text)) if node_text[i] in pathway_set]
        adjust_text(texts,arrowprops=dict(arrowstyle='-',color='#CFCFD0'))
        plt.savefig(os.path.join(outdir,'GO_plot_static.pdf'),bbox_inches='tight');plt.close()

            

        







def mutation_analysis(mode,burden,mutation,output,n_sample_cutoff=10,gene_column='gene',genes_to_plot=None):
    '''
    Run mutation association analysis with neoantigen burden

    :param mode: string, either 'compute' or 'plot', compute will compute all the association, plot will visualize certain assciation in box plot
    :param burden: pandas dataframe, the burden result file generated by T pipeline, and sample ID should be consistent with the sample ID in mutation dataframe index
    :param mutation: pandas dataframe, the mutation data, make sure sample ID column is the index
    :param output: string, the path in which the result will be generated, if compute mode, a txt file will be generated, if plot mode, a plot will be generated
    :param n_sample_cutoff: int, default is 10, if a mutation occur in < n_sample_cutoff samples, it won't be considered
    :param gene_column: string, default is gene, in the mutation df, which column represents the gene
    :param gene_to_plot: list, each item is a gene name, for the plot mode

    Example::

        # compute mode
        snaf.mutation_analysis(mode='compute',burden=burden,mutation=mutation,output='result/stage3_mutation.txt')
        # plot mode
        snaf.mutation_analysis(mode='plot',burden=burden,mutation=mutation,output='result/stage3_mutation_CAMKK2.pdf',genes_to_plot=['CAMKK2'])
    '''
    if mode == 'compute':
        burden = burden.loc[:,burden.columns.isin(mutation.index)].iloc[-1]
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
            burden_df = burden.loc[:,burden.columns.isin(mutation.index)].iloc[-1].to_frame()
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
    '''
    Run survival analysis based on neoantigen burden, the `burden`and `survival`must be pandas dataframe, and make sure the 
    sample ID are consistent with each other, we only consider samples in `burden` that have clinical survival data.

    :param burden: pandas dataframe, the burden result file generated by T pipeline
    :param survival: pandas dataframe, the survival data
    :param n: int, how many groups to stratigy, support 2 (median), 3 (33%,66%), 4 (quantile)
    :param stratification_plot: string, the path and name for the stratification plot
    :param survival_plot: string, the path plus name for the survival plot
    :param survival_duration: string, which column in survival data represent the duration
    :param survival_event: string, which column in survival data represent the event

    Example::

        snaf.survival_analysis(burden,survival,n=2,stratification_plot='result/stage2_stratify.pdf',survival_plot='result/stage2_survival.pdf')

    '''

    burden = burden.loc[:,burden.columns.isin(survival.index)].iloc[-1]
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
        quantiles = burden.quantile([0.25,0.5,0.75]).values
        iqr = quantiles[2] - quantiles[0]
        upper_bound = quantiles[2] + 1.5*iqr
        lower_bound = quantiles[0] - 1.5*iqr
        identity_col = []
        for item in burden:
            if item > upper_bound:
                identity_col.append('outlier')
            elif item > quantiles[1] and item <= upper_bound:
                identity_col.append('high')
            elif item >= lower_bound and item <= quantiles[1]:
                identity_col.append('low')
            else:
                identity_col.append('outlier')
    burden_output = burden.to_frame()
    burden_output['identity'] = identity_col
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
    ax.legend(current_handles,new_labels,bbox_to_anchor=(1,1),loc='upper left',frameon=False)
    results = logrank_test(low_os[survival_duration],high_os[survival_duration],low_os[survival_event],high_os[survival_event])
    ax.text(x=1000,y=0.05,s='Log-rank test: p-value is {:.2f}'.format(results.p_value),weight='bold')
    plt.savefig(survival_plot,bbox_inches='tight');plt.close()
    return burden_output,quantiles


def report_candidates(jcmq,df,sample,outdir,remove_quote=True,metrics={'netMHCpan_el':'binding_affinity','deepimmuno_immunogenicity':'immunogenicity'},
                      criterion=[('netMHCpan_el',0,'<=',2),('deepimmuno_immunogenicity',1,'==','True'),]):
    '''
    this function will report a txt file for all the T antigen candidate for a specific tumor sample.

    :param jcmq: the JunctionCountMatrixQuery object, usually you need to deserialize after_prediction.p file produced by T antigen pipeline
    :param df: DataFrame,the frequency table generated by T antigen pipeline
    :param sample: String, the sample name to report T antigen candiate
    :param outdir: string, the path to the output directory
    :param remove_quote: boolean, whether to remove the quotation or not, as one column in frequency table df is list, when loaded in memory using pandas, it will be added a quote, we can remove it
    :param metrics: dict, the key is the metric we registered to each peptide-hla, values is the corresponding name we want to output in the report for each metric
    :param criterion: list with each element as a tuple, this is a very specialized expression to specificy the criterion we are going to use to only report peptide-hla that meet the criterion.
                      each tuple has four element, the first element is the metric name we registered, the second is either 0 or 1, 0 means looking for the continuous value score,
                      1 means looking for the string identity value. See below note for further explanation. the third is an operator, and the fourth is the cutoff.

    .. note::

        By default, we run netMHCpan or MHCflurry for MHC binding prediction, and DeepImmuno for immunogenicity prediction. When we add the predicted value for
        each peptide-hla pair, we format the prediction result into a dataframe compliant with something like below:

          peptide     mer     hla       score    identity
         AAAAAAAAA     9   HLA-A*01:01   0.3       SB  

        here is the dataframe for MHC binding prediction, we have a column recording continous value (score), and another column (identity) recording whether it is
        a SB (Strong Binder) or WB (Weaker Binder). Those are all useful information coming out of the predictor. Now we just need to register this df to the NeoJunction
        class by running ``nj.enhanced_peptide.register_attr(df,'netMHCpan_el)``. We can keep adding scores to the NeoJunction class. This explain what the 0/1 in criterion
        mean.

    Example::

        snaf.report_candidates(jcmq,df,'SRR067783.bed','result',True)
        # a txt file will be written to the outdir you specified

    '''
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    df = df.copy()
    from ast import literal_eval
    if remove_quote:
        df['samples'] = [literal_eval(item) for item in df['samples']]
    # build dict for retrieving specificity scores and occurence
    df_score = df.filter(like='tumor_specificity',axis=1)
    score_dict = {}   # {tumor_specificity_mle:{aa1_uid1:0.5,aa2_uid2:0.6},tumor_specifity_bayesian:{}...}
    for col in df_score.columns.tolist() + ['n_sample','coord','symbol']:
        tmp = df[col].to_dict()
        score_dict[col] = tmp
    # start to process jcmq
    junction_count_dict = jcmq.junction_count_matrix[sample].to_dict()
    col_index = jcmq.subset.columns.tolist().index(sample)
    results = jcmq.results[0]
    hlas = jcmq.results[1]
    selected_hla = hlas[col_index]
    with open(os.path.join(outdir,'T_antigen_candidates_{}.txt'.format(sample)),'w') as f:
        f.write('sample\tpeptide\tuid\tjunction_count\tphase\tevidences\thla\t')
        metrics_stream = '\t'.join(list(metrics.values())) + '\t'
        f.write(metrics_stream)
        score_stream = '\t'.join(list(score_dict.keys())) + '\n'
        f.write(score_stream)
        for item,samples in zip(df.index,df['samples']):
            if sample in samples:
                stream = '{}\t'.format(sample)
                aa,uid = item.split(',')
                jc = junction_count_dict[uid]
                stream += '{}\t{}\t{}\t'.format(aa,uid,jc)
                row_index = jcmq.subset.index.tolist().index(uid)
                nj = deepcopy(results[row_index])
                nj.enhanced_peptides = nj.enhanced_peptides.filter_based_on_hla(selected_hla=selected_hla)
                ep = nj.enhanced_peptides.filter_based_on_criterion(criterion,True)  # only report valid hla
                origin = ep[len(aa)][aa]['origin']
                stream += '{}\t{}\t'.format(origin[2],origin[3])  # phase, evidences
                for hla in ep[len(aa)][aa].keys():
                    if hla != 'origin':
                        stream += '{}\t'.format(hla)
                        for k,v in metrics.items():
                            s = nj.enhanced_peptides[len(aa)][aa][hla][k][0]
                            stream += '{}\t'.format(s)
                        for k,v in score_dict.items():
                            s = v[item]
                            stream += '{}\t'.format(s)
                        f.write(stream.rstrip('\t') + '\n')
                        stream = '\t'.join(stream.split('\t')[:6]) + '\t'



    


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


def add_gene_symbol_frequency_table(df,remove_quote=True):
    '''
    This function will convert the ENSG id to gene sysmbol and add a column in place for your dataframe

    :param df: the input df, make sure the index is the uid
    :param remove_quote: bool, depending on how your df is loaded, if directly read from the disk, certain column will contain quotation, whether to remove it or not

    :return df: the output df, with added gene symbol column

    Example::

        add_gene_symbol_frequency_table(df=df,remove_quote=False)
    '''
    # the index has to be comma separated string
    from ast import literal_eval
    if remove_quote:
        df['samples'] = [literal_eval(item) for item in df['samples']]
    ensg_list = [item.split(',')[1].split(':')[0] for item in df.index] 
    symbol_list = ensemblgene_to_symbol(ensg_list,'human') 
    df['symbol'] = symbol_list
    return df


def reformat_frequency_table(df,remove_quote=True):
    # the index has to be comma separated string
    from ast import literal_eval
    if remove_quote:
        df['samples'] = [literal_eval(item) for item in df['samples']]
    sequence_io = []
    for row in df.itertuples():
        for item in row.samples:
            sequence_io.append((row.Index,item,1))
    result = pd.DataFrame.from_records(sequence_io,columns=['id','sample','value'])
    result = result.groupby(by=['id','sample'])['value'].mean().unstack(fill_value=0)
    return result
        

def ensemblgene_to_symbol(query,species):
    '''
    Examples::
        from sctriangulate.preprocessing import GeneConvert
        converted_list = GeneConvert.ensemblgene_to_symbol(['ENSG00000010404','ENSG00000010505'],species='human')
    '''
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


def analyze_neoantigens(freq_path,junction_path,total_samples,outdir,fasta=False,mers=None,columns=None,cutoffs=(0.1,0.9),junction_bl=0.1):
    '''
    Users need to execute this function to prepare input for the t antigen dash viewer

    :param freq_path: string, the path to the frequency file in the result folder, with uid appended
    :param junction_path: string, the path to the burden0 file in the result folder
    :param total_samples: int, the total number of samples, will be used for determine shared versus unique neoantigen
    :param outdir: string, where the output will go into
    :param fasta: bool, whether to output the fasta file as well to check the motif 
    :param mers: None or list, if [9,10], it will generate mer9 and mer10 for fasta and separate dataframe.
    :param columns: None or string, which column to include in the output df
    :param cutoffs: tuple, the cutoffs to determine shared versus unique neoantigen
    :param junction_bl: float, the cutff for parent junction, if junction occur < junction_bl sample, we don't consider neoantigens from those junction

    Example::

        # just output for viewer
        snaf.downstream.analyze_neoantigens(freq_path='result/frequency_stage2_verbosity1_uid.txt',junction_path='result/burden_stage0.txt',total_samples=472,outdir='result',mers=None,fasta=False)
        # to run MEME for motif after that
        snaf.downstream.analyze_neoantigens(freq_path='result/frequency_stage2_verbosity1_uid.txt',junction_path='result/burden_stage0.txt',total_samples=472,outdir='result',mers=[9,10],fasta=True)
    '''
    freq = pd.read_csv(freq_path,sep='\t',index_col=0)
    junction = pd.read_csv(junction_path,sep='\t',index_col=0)['mean'].to_dict()
    freq['uid'] = [item.split(',')[-1] for item in freq.index]
    freq.index = [item.split(',')[0] for item in freq.index]
    freq['mean_percent_samples_junction_present'] = freq['uid'].map(junction).values
    freq['actual_percent_samples_neoantigen_present'] = freq['n_sample'] / total_samples
    freq.drop(columns='samples',inplace=True)
    identity = []
    for e,o in zip(freq['mean_percent_samples_junction_present'],freq['actual_percent_samples_neoantigen_present']):
        if o < cutoffs[0] * e:
            identity.append('low')
        elif o > cutoffs[1] * e:
            identity.append('high')
        else:
            identity.append('medium')
    freq['identity'] = identity
    freq = freq.loc[np.logical_not(freq.index.duplicated()),:]
    freq = freq.loc[freq['identity']!='medium',:]
    freq = freq.loc[freq['mean_percent_samples_junction_present']>junction_bl,:]
    freq['length'] = [len(item) for item in freq.index]
    if columns is None:
        columns = []
    selected_columns = ['n_sample','uid','mean_percent_samples_junction_present','actual_percent_samples_neoantigen_present','identity','length'] + columns
    freq = freq.loc[:,selected_columns]
    sns.regplot(data=freq,x='mean_percent_samples_junction_present',y='actual_percent_samples_neoantigen_present',fit_reg=False,x_jitter=0.01,y_jitter=0.01,scatter_kws={'s':5})
    plt.savefig(os.path.join(outdir,'shared_vs_unique_neoantigen_all.pdf'),bbox_inches='tight')
    plt.close()
    if mers is None:
        freq.to_csv(os.path.join(outdir,'shared_vs_unique_neoantigen_all.txt'),sep='\t')
    else:
        for mer in mers:
            freq_mer = freq.loc[freq['length']==mer,:]
            freq_mer.to_csv(os.path.join(outdir,'shared_vs_unique_neoantigen_mer{}.txt'.format(mer)),sep='\t')
            if fasta:
                with open(os.path.join(outdir,'mer{}_high.fasta'.format(mer)),'w') as f1, open(os.path.join(outdir,'mer{}_low.fasta'.format(mer)),'w') as f2:  
                    for identity,sub_df in freq_mer.groupby(by='identity'):
                        if identity == 'high':
                            for item in sub_df.index:
                                f1.write('>{}\n{}\n'.format(item,item))
                        else:
                            for item in sub_df.index:
                                f2.write('>{}\n{}\n'.format(item,item))  
                     









    