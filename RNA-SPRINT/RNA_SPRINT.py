#!/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/decoupler_env/bin/python3.8

import decoupler as dc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os,sys
from tqdm import tqdm
from adjustText import adjust_text
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# you only need to change this captical variable, a file called mdt_estimate.tsv will show estimated RBP activity for each sample 
EVENT_ANNOTATION_FILE_PATH = '../TCGA_melanoma/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt'
PROPORTION_CUTOFF = 0.75
VARIATION_CUTOFF = 0.1
K562_PRIOR_NETWORK_PATH = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/HepG2_inference/prior.tsv'
OUTPUT_DIR = '.'

# preprocessing to get necessary input files, code derived from https://github.com/frankligy/splice-inferelator/blob/main/scripts/prepare_inference.py
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

def median_impute(x):
    med = np.ma.median(np.ma.masked_invalid(x.values))
    result = np.nan_to_num(x.values,nan=med)
    return result

ee = pd.read_csv(EVENT_ANNOTATION_FILE_PATH,sep='\t',index_col='UID').iloc[:,10:]
ee.index = [disentangle_uid(item,add_symbol=True) for item in ee.index]
prop = 1- np.count_nonzero(ee.isna().values,axis=1)/ee.shape[1]  # 101974 rows x 472 columns
ee = ee.loc[prop>PROPORTION_CUTOFF,:]   # 59134 rows x 472 columns
original_columns = ee.columns
ee = ee.apply(func=median_impute,axis=1,result_type='expand')
ee.columns = original_columns
cv = ee.apply(func=lambda x:x.values.std()/x.values.mean(),axis=1,result_type='reduce').values
ee = ee.loc[cv>VARIATION_CUTOFF,:]
ee.columns = original_columns   # 34683 rows x 472 columns
ee.T.to_csv('expr.tsv',sep='\t')


expr_path = 'expr.tsv'
prior_path = K562_PRIOR_NETWORK_PATH
root_path = OUTPUT_DIR
expr = pd.read_csv(expr_path,sep='\t',index_col=0)  # 472 rows x 34683 columns
prior = pd.read_csv(prior_path,sep='\t',index_col=0)  # 73157 rows x 221 columns
common_events = list(set(expr.columns).intersection(set(prior.index)))
expr = expr.loc[:,common_events]   # 472 rows x 17836 columns
prior = prior.loc[common_events,:]  # 17836 rows x 221 columns
net = prior.stack().reset_index()
net.columns = ['target','source','weight']
net = net.loc[net['weight']!=0,:]

estimate = dc.run_mdt(mat=expr,net=net)
estimate.to_csv(os.path.join(root_path,'mdt_estimate.txt'),sep='\t')