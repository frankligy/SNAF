#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os,subprocess,sys
from scipy.stats import norm
from tqdm import tqdm
import pickle
import plotly.graph_objects as go
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import os,subprocess,sys
from scipy.stats import norm
from tqdm import tqdm
import pickle
import plotly.graph_objects as go
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt

def find_last_exon(exon):
    if exon.startswith('E'): # it is an exon
        try:
            assert int(exon[1:])-1 >= 1
        except AssertionError:   # means this exon, E1, is already the first one
            last = exon
        else:
            last = 'I' + str(int(exon[1:])-1)
    elif exon.startswith('I'): # it is an intron
        last = 'E' + exon[1:]
    return last

def find_next_exon(exon):
    '''
    One caveat for this function, E5 might be the last exon in the transcript
    '''
    if exon.startswith('E'):  # it is an exon
        next_ = 'I' + exon[1:]
    elif exon.startswith('I'):  # it is an intron
        next_ = 'E' + str(int(exon[1:]) + 1)
    return next_

def return_coords(big_dict,gene,query,start,end):
    '''
    big_dict stores the coord for every exon, not subexon
    gene is the ENSG id
    query is the query exon
    start is the exon that precedes the query exon, getting from find_last_exon function
    end is the exon that follows the query exon, getting from find_next_exon function
    '''
    global error   # a exon that is not present in db, 
    strand = big_dict[gene][0]
    chromosome = big_dict[gene][1]
    if strand == '+':
        try:
            final_1 = big_dict[gene][2][start][0]    # left-most
        except:
            error += 1
            return [1]
        try:
            final_2 = big_dict[gene][2][end][1]   # right-most
        except KeyError:   # E5 is the last exon, so if I5, it will throw a KeyError
            try:
                final_2 = big_dict[gene][2][query][1]   
            except:
                error += 1
                return [1]
    else:
        try:
            final_1 = big_dict[gene][2][end][0]
        except KeyError:
            try:
                final_1 = big_dict[gene][2][query][0]
            except:
                error += 1
                return [1]
        try:
            final_2 = big_dict[gene][2][start][1]
        except:
            error += 1
            return [1]
    '''return final1 always smaller than final2, because we are talking about the absolute position in positive strand'''
    return chromosome,strand,final_1,final_2

def disentangle_uid(uid,add_stable_id=False,add_symbol=False):
    fg,bg = uid.split('|')
    gene_symbol = fg.split(':')[0]
    stable_id = fg.split(':')[1]
    exons = ':'.join(fg.split(':')[2:])
    output = exons
    if add_stable_id:
        output = stable_id + ':' + output
    if add_symbol:
        output = gene_symbol + ':' + output
    return output

if __name__ == '__main__':

    sf = 'PCBP2'
    uid = 'PCK2:E10.1-E11.1'

    group_file_combine_K562 = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/das_combine/groups.{}_shRNA_K562.txt'.format(sf)
    group_file_combine_HepG2 = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/das_combine/groups.{}_shRNA_HepG2.txt'.format(sf)
    group_file_within_K562 = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/das_within/groups.{}_shRNA_K562.txt'.format(sf)
    group_file_within_HepG2 = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/das_within/groups.{}_shRNA_HepG2.txt'.format(sf)
    das_file_combine_K562 = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_K562/{}.txt'.format(sf)
    das_file_combine_HepG2 = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_HepG2/{}.txt'.format(sf)
    psi_matrix_file_combine = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/limma_bcAltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt'
    psi_matrix_file_within = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt'

    '''first part, strip plot'''

    fig,axes = plt.subplots(nrows=1,ncols=2,gridspec_kw={'wspace':0.5})
    # draw axes[0], K562
    tmp = pd.read_csv(group_file_combine_K562,sep='\t',header=None,index_col=0)
    exp_columns = tmp.loc[tmp[2]=='exp'].index.tolist()
    control_columns = tmp.loc[tmp[2]=='control'].index.tolist()
    tmp = pd.read_csv(group_file_within_K562,sep='\t',header=None,index_col=0)
    match_columns = tmp.loc[tmp[2]=='control'].index.tolist()

    tmp = pd.read_csv(psi_matrix_file_combine,sep='\t',index_col='UID').iloc[:,10:]
    tmp['new_uid'] =[disentangle_uid(item,add_symbol=True) for item in tmp.index]
    exp_values = tmp.loc[tmp['new_uid']==uid,exp_columns].values.squeeze()
    control_values = tmp.loc[tmp['new_uid']==uid,control_columns].values.squeeze()
    df = pd.concat(objs=[pd.DataFrame({'identity':['exp']*len(exp_values),'value':exp_values}),pd.DataFrame({'identity':['control']*len(control_values),'value':control_values})],axis=0)
    sns.stripplot(x='identity',y='value',data=df,ax=axes[0])


    tmp = pd.read_csv(das_file_combine_K562,sep='\t',index_col=0)
    tmp['new_uid'] = [disentangle_uid(item,add_symbol=True) for item in tmp.index]
    metrics = ','.join([str(round(item,4)) for item in tmp.loc[tmp['new_uid']==uid,['dPSI','rawp','adjp']].values.squeeze()])
    axes[0].text(x=0.5,y=1,s=metrics,fontsize=4,va='top',ha='center',transform=axes[0].transAxes)

    # draw axes[1], HepG2
    tmp = pd.read_csv(group_file_combine_HepG2,sep='\t',header=None,index_col=0)
    exp_columns = tmp.loc[tmp[2]=='exp'].index.tolist()
    control_columns = tmp.loc[tmp[2]=='control'].index.tolist()
    tmp = pd.read_csv(group_file_within_HepG2,sep='\t',header=None,index_col=0)
    match_columns = tmp.loc[tmp[2]=='control'].index.tolist()

    tmp = pd.read_csv(psi_matrix_file_combine,sep='\t',index_col='UID').iloc[:,10:]
    tmp['new_uid'] =[disentangle_uid(item,add_symbol=True) for item in tmp.index]
    exp_values = tmp.loc[tmp['new_uid']==uid,exp_columns].values.squeeze()
    control_values = tmp.loc[tmp['new_uid']==uid,control_columns].values.squeeze()
    df = pd.concat(objs=[pd.DataFrame({'identity':['exp']*len(exp_values),'value':exp_values}),pd.DataFrame({'identity':['control']*len(control_values),'value':control_values})],axis=0)
    sns.stripplot(x='identity',y='value',data=df,ax=axes[1])

    tmp = pd.read_csv(das_file_combine_HepG2,sep='\t',index_col=0)
    tmp['new_uid'] = [disentangle_uid(item,add_symbol=True) for item in tmp.index]
    metrics = ','.join([str(round(item,4)) for item in tmp.loc[tmp['new_uid']==uid,['dPSI','rawp','adjp']].values.squeeze()])
    axes[1].text(x=0.5,y=1,s=metrics,fontsize=4,va='top',ha='center',transform=axes[1].transAxes)

    # print
    fig.suptitle('{}->{} functional left-K562 right-HepG2'.format(sf,uid))
    plt.savefig('check.pdf',bbox_inches='tight');plt.close()

    '''second, interval of the splicing event'''

    ee_path = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt'
    event_and_alt = pd.read_csv(ee_path,sep='\t').loc[:,['UID','AltExons']].set_index('UID').squeeze()  # series
    event_and_alt.index = [disentangle_uid(item,add_symbol=True) for item in event_and_alt.index]
    event_and_alt = event_and_alt.loc[np.logical_not(event_and_alt.index.duplicated())]
    reg_exon_list = event_and_alt.to_dict()[uid].split('|')

    with open('../big_dict.p','rb') as f:
        big_dict = pickle.load(f)

    for exon in reg_exon_list:
        gene = exon.split(':')[0]
        query = exon.split(':')[1].split('.')[0]
        start = find_last_exon(query)
        end = find_next_exon(query)
        chrom, strand, final_1, final_2 = return_coords(big_dict,gene,query,start,end)
        print('{} interval:\n{}:{}-{}'.format(exon,chrom,final_1,final_2))

