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

def mutation_analysis(mode,burden,mutation,output,n_sample_cutoff=10,gene_column='gene',genes_to_plot=None):
    '''
    Run mutation association analysis with neoantigen burden

    :param mode: string, either 'compute' or 'plot', compute will compute all the association, plot will visualize certain assciation in box plot
    :param burden: pandas dataframe, the burden result file generated by T pipeline
    :param mutation: pandas dataframe, the mutation data
    :param n_sample_cutoff: int, default is 10, if a mutation occur in < n_sample_cutoff samples, it won't be considered
    :param gene_column: string, default is gene, in the mutation df, which column represents the gene
    :param gene_to_plot: list, each item is a gene name, for the plot mode

    Example::

        # compute mode
        snaf.mutation_analysis(mode='compute',burden=burden,mutation=mutation,output='result/stage3_mutation.txt')
        # plot mode
        snaf.mutation_analysis(mode='plot',burden=burden,mutation=mutation,output='result/stage3_mutation_CAMKK2.pdf',genes_to_plot=['CAMKK2'])
    '''
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
    '''
    Run survival analysis based on neoantigen burden

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
    from ast import literal_eval
    if remove_quote:
        df['samples'] = [literal_eval(item) for item in df['samples']]
    # build dict for retrieving specificity scores and occurence
    df_score = df.filter(like='tumor_specificity',axis=1)
    score_dict = {}   # {tumor_specificity_mle:{aa1_uid1:0.5,aa2_uid2:0.6},tumor_specifity_bayesian:{}...}
    for col in df_score.columns.tolist() + ['n_sample']:
        tmp = df[col].to_dict()
        score_dict[col] = tmp
    # start to process jcmq
    col_index = jcmq.subset.columns.tolist().index(sample)
    results = jcmq.results[0]
    hlas = jcmq.results[1]
    selected_hla = hlas[col_index]
    with open(os.path.join(outdir,'T_antigen_candidates_{}.txt'.format(sample)),'w') as f:
        f.write('sample\tpeptide\tuid\thla\t')
        metrics_stream = '\t'.join(list(metrics.values())) + '\t'
        f.write(metrics_stream)
        score_stream = '\t'.join(list(score_dict.keys())) + '\n'
        f.write(score_stream)
        for item,samples in tqdm(zip(df.index,df['samples']),total=df.shape[0]):
            if sample in samples:
                stream = '{}\t'.format(sample)
                aa,uid = item.split(',')
                stream += '{}\t{}\t'.format(aa,uid)
                row_index = jcmq.subset.index.tolist().index(uid)
                nj = deepcopy(results[row_index])
                nj.enhanced_peptides = nj.enhanced_peptides.filter_based_on_hla(selected_hla=selected_hla)
                ep = nj.enhanced_peptides.filter_based_on_criterion(criterion,True)  # only report valid hla
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
                        stream = '\t'.join(stream.split('\t')[:3]) + '\t'



    


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
    result = out['out']['symbol'].fillna('unknown_gene').tolist()
    try:
        assert len(query) == len(result)
    except AssertionError:    # have duplicate results
        df = out['out']
        df_unique = df.loc[~df.index.duplicated(),:]
        result = df_unique['symbol'].fillna('unknown_gene').tolist()
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
                     









    