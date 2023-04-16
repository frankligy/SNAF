#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os,subprocess,sys
from scipy.stats import norm,rankdata
from tqdm import tqdm
import pickle


def weighted_z_score_combine(pval_list,weights=None):
    zscore_list = np.array([norm.ppf(pval,0,1) for pval in pval_list])
    if weights is None:
        weights = np.full(shape=len(zscore_list),fill_value=1/len(zscore_list))
    combined_zscore = np.inner(weights,zscore_list)/np.sqrt(np.sum(weights))
    combined_pval = norm.cdf(combined_zscore,0,1)
    return combined_pval
    

def percentile_combine(sf,series_list,cutoff=0.05):
    event_union = set()
    sorted_frames_list = []  
    for series in series_list:
        frames = series.loc[series<cutoff].sort_values().to_frame(name='rawp').reset_index()
        sorted_frames_list.append(frames)
        event_union = event_union.union(set(frames['index'].values))
    results_list = []
    for event in event_union:
        percentile_list = []
        for frames in sorted_frames_list:
            try:
                percentile = frames.loc[frames['index']==event,:].index.values[0]/frames.shape[0]
            except IndexError:  # empty Index, empty 1d array, meaning it is not significant
                continue
            percentile_list.append(percentile)
        if len(percentile_list) == 0:
            combined_percentile = 1
        else:
            combined_percentile = np.array(percentile_list).mean()
        results_list.append((sf,event,1-combined_percentile))   # higher, the better
    return results_list

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

def fisher_method(p1,p2):
    from scipy.stats import chi2
    k = 2
    chi_square_statistics = -2 * (np.log(p1) + np.log(p2))
    pval = chi2.sf(chi_square_statistics,df=2*k)
    return pval


# function evidence rank combine
os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/shRNA_K562')
shRNA_K562_sfs = subprocess.run('for file in *.txt; do echo $(basename $file .txt); done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/shRNA_HepG2')
shRNA_HepG2_sfs = subprocess.run('for file in *.txt; do echo $(basename $file .txt); done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]    
os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/CRISPR_K562')
CRISPR_K562_sfs = subprocess.run('for file in *.txt; do echo $(basename $file .txt); done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/CRISPR_HepG2')
CRISPR_HepG2_sfs = subprocess.run('for file in *.txt; do echo $(basename $file .txt); done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]    
all_sfs = list(set(shRNA_K562_sfs).union(set(shRNA_HepG2_sfs)).union(set(CRISPR_K562_sfs)).union(set(CRISPR_HepG2_sfs))) # 299
d = 0.1
p = 0.05


function_prior = []
for sf in tqdm(shRNA_HepG2_sfs):
    # combine
    os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_HepG2')
    cwd = os.getcwd()
    psi = pd.read_csv(cwd + '/' + sf + '.txt',sep='\t',index_col=0)
    psi.index = [disentangle_uid(item,add_symbol=True) for item in psi.index]
    psi = psi.loc[np.logical_not(psi.index.duplicated()),:]
    psi = psi.loc[(psi['dPSI'].abs()>d)&(psi['rawp']<p),:]
    psi['rank_d'] = psi.shape[0] - rankdata(psi['dPSI'].abs().values)
    psi['rank_p'] = rankdata(psi['rawp'].values)
    psi['rank'] = (psi['rank_d'] + psi['rank_p'])/2
    frame_combine = psi['rank'].sort_values().to_frame(name='rank').reset_index()
    frame_combine['percentile'] = frame_combine.index.values / frame_combine.shape[0]
    # within
    os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/shRNA_HepG2')
    cwd = os.getcwd()
    psi = pd.read_csv(cwd + '/' + sf + '.txt',sep='\t',index_col=0)
    psi.index = [disentangle_uid(item,add_symbol=True) for item in psi.index]
    psi = psi.loc[np.logical_not(psi.index.duplicated()),:]
    psi = psi.loc[(psi['dPSI'].abs()>d)&(psi['rawp']<p),:]
    psi['rank_d'] = psi.shape[0] - rankdata(psi['dPSI'].abs().values)
    psi['rank_p'] = rankdata(psi['rawp'].values)
    psi['rank'] = (psi['rank_d'] + psi['rank_p'])/2
    frame_within = psi['rank'].sort_values().to_frame(name='rank').reset_index()
    frame_within['percentile'] = frame_within.index.values / frame_within.shape[0]
    # inner join
    frame = frame_combine.join(frame_within.set_index('index'),on='index',how='inner',lsuffix='_combine',rsuffix='_within')
    frame['meta_percentile'] = (frame['percentile_combine'] + frame['percentile_within']) / 2
    frame = frame.sort_values(by='meta_percentile')
    frame = frame.set_index(keys=np.arange(frame.shape[0]))
    for row in frame.itertuples():
        percentile = row.Index / frame.shape[0]
        edge = (sf,row.index,1-percentile)
        function_prior.append(edge)
    





# function_prior = []
# for sf in tqdm(all_sfs):
#     series_list = []
#     # go to shRNA_K562
#     os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/shRNA_K562')
#     cwd = os.getcwd()
#     try:
#         psi = pd.read_csv(cwd + '/' + sf + '.txt',sep='\t',index_col=0)
#         psi.index = [':'.join(item.split('|')[0].split(':')[1:]) for item in psi.index]
#         psi = psi.loc[np.logical_not(psi.index.duplicated()),:]
#     except FileNotFoundError:
#         pass
#     else:
#         series_list.append(psi['rawp'])
#     # go to shRNA_HepG2
#     os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/shRNA_HepG2')
#     cwd = os.getcwd()
#     try:
#         psi = pd.read_csv(cwd + '/' + sf + '.txt',sep='\t',index_col=0)
#         psi.index = [':'.join(item.split('|')[0].split(':')[1:]) for item in psi.index]
#         psi = psi.loc[np.logical_not(psi.index.duplicated()),:]
#     except FileNotFoundError:
#         pass
#     else:
#         series_list.append(psi['rawp'])
#     # go to CRISPR_K562
#     os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/CRISPR_K562')
#     cwd = os.getcwd()
#     try:
#         psi = pd.read_csv(cwd + '/' + sf + '.txt',sep='\t',index_col=0)
#         psi.index = [':'.join(item.split('|')[0].split(':')[1:]) for item in psi.index]
#         psi = psi.loc[np.logical_not(psi.index.duplicated()),:]
#     except FileNotFoundError:
#         pass
#     else:
#         series_list.append(psi['rawp'])    
#     # go to CRISPR_HepG2
#     os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/CRISPR_HepG2')
#     cwd = os.getcwd()
#     try:
#         psi = pd.read_csv(cwd + '/' + sf + '.txt',sep='\t',index_col=0)
#         psi.index = [':'.join(item.split('|')[0].split(':')[1:]) for item in psi.index]
#         psi = psi.loc[np.logical_not(psi.index.duplicated()),:]
#     except FileNotFoundError:
#         pass
#     else:
#         series_list.append(psi['rawp'])    
#     # run percentile combine
#     results_list = percentile_combine(sf,series_list)
#     function_prior.extend(results_list)
# with open('/data/salomonis-archive/FASTQ/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/prior/function_prior.p','wb') as f:
    # pickle.dump(function_prior,f,protocol=pickle.HIGHEST_PROTOCOL)

# with open('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/prior/function_prior.p','rb') as f:
#     function_prior = pickle.load(f)
prior_function = pd.DataFrame()
tmp = list(zip(*function_prior))
prior_function['sf'] = tmp[0]
prior_function['event'] = tmp[1]
prior_function['value'] = tmp[2]
prior_function = prior_function.groupby(by=['event','sf'])['value'].mean().unstack(fill_value=0)
print(prior_function)
prior_function.to_csv('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/prior/functional_prior/shRNA_HepG2_consensus.txt',sep='\t')

















