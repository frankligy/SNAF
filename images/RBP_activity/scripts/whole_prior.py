#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os,subprocess,sys
from scipy.stats import norm
from tqdm import tqdm
import pickle

func_prior = pd.read_csv('../altanalyze_output/prior/functional_prior/shRNA_HepG2_consensus.txt',sep='\t',index_col=0).to_dict(orient='dict')
bind_prior = pd.read_csv('../altanalyze_output/prior/binding_prior/HepG2_prior.txt',sep='\t',index_col=0).to_dict(orient='dict')

func_prior_dict = {}
for col,row in func_prior.items():
    for index,entry in row.items():
        func_prior_dict[(col,index)] = entry
bind_prior_dict = {}
for col,row in bind_prior.items():
    for index,entry in row.items():
        bind_prior_dict[(col,index)] = entry


def merge_by_rule(func_prior_dict,bind_prior_dict,key):
    both = 1
    func = 0.49
    bind = 0.51
    none = 0
    func_evidence = func_prior_dict.get(key,0)
    bind_evidence = bind_prior_dict.get(key,0)
    if func_evidence > 0 and bind_evidence > 0:
        return both
    elif func_evidence > 0 and bind_evidence == 0:
        return func
    elif func_evidence ==0 and bind_evidence > 0:
        return bind
    else:
        return none


# merge the two dict
whole_prior_dict = {key:merge_by_rule(func_prior_dict,bind_prior_dict,key) for key in set(func_prior_dict).union(set(bind_prior_dict))}
whole_prior_list = []
for key,value in whole_prior_dict.items():
    whole_prior_list.append((key[0],key[1],value))

# to_df
prior_df = pd.DataFrame()
tmp = list(zip(*whole_prior_list))
prior_df['sf'] = tmp[0]
prior_df['event'] = tmp[1]
prior_df['value'] = tmp[2]
prior_df = prior_df.groupby(by=['event','sf'])['value'].mean().unstack(fill_value=0)
print(prior_df)
prior_df.to_csv('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/prior/whole_prior/shRNA_HepG2_rule.txt',sep='\t')

















