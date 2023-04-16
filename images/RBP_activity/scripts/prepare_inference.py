#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os,subprocess,sys
from scipy.stats import norm
from tqdm import tqdm
import pickle

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

# if using split eventannotation file, the additional column have alrady been removed, otherwise, make sure to preprocess
ee = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/AML/Leucegene/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t',index_col='UID').iloc[:,10:]
ee.index = [disentangle_uid(item,add_symbol=True) for item in ee.index]
prop = 1- np.count_nonzero(ee.isna().values,axis=1)/ee.shape[1]
print(ee)
ee = ee.loc[prop>0.75,:]  
print(ee)
original_columns = ee.columns
ee = ee.apply(func=median_impute,axis=1,result_type='expand')
ee.columns = original_columns
cv = ee.apply(func=lambda x:x.values.std()/x.values.mean(),axis=1,result_type='reduce').values
ee = ee.loc[cv>0.1,:]
ee.columns = original_columns   # [35573 rows x 446 columns]
print(ee)

# no split, just run as a whole
'''
expr: sample * event
prior: event * sf
sf_name: sf
'''
expr = ee.T   
prior = pd.read_csv('../altanalyze_output/prior/whole_prior/shRNA_K562_rule.txt',sep='\t',index_col=0)   # 28733 commonly present events
expr.to_csv('../altanalyze_output/leucegene_inference/expr.tsv',sep='\t')
prior.to_csv('../altanalyze_output/leucegene_inference/prior.tsv',sep='\t')
prior.columns.to_series().to_csv('../altanalyze_output/leucegene_inference/sf_name.tsv',sep='\t',index=None,header=None)























