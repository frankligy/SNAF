#!/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/RNA-SPRINT/rna_sprint_env/bin/python3.8

import decoupler as dc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os,sys
from tqdm import tqdm
import matplotlib as mpl
import argparse

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


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

def main(args):
    EVENT_ANNOTATION_FILE_PATH = args.splicing
    PROPORTION_CUTOFF = args.prop
    VARIATION_CUTOFF = args.cv
    K562_PRIOR_NETWORK_PATH = args.prior
    OUTPUT_DIR = args.outdir
    meta = args.meta

    if not os.path.exists(OUTPUT_DIR):
        os.mkdir(OUTPUT_DIR)
    root_path = OUTPUT_DIR

    ee = pd.read_csv(EVENT_ANNOTATION_FILE_PATH,sep='\t',index_col='UID').iloc[:,10:]
    print('loaded splicing matrix: {}'.format(ee.shape))
    ee.index = [disentangle_uid(item,add_symbol=True) for item in ee.index]
    prop = 1- np.count_nonzero(ee.isna().values,axis=1)/ee.shape[1]  # 101974 rows x 472 columns
    ee = ee.loc[prop>PROPORTION_CUTOFF,:]   # 59134 rows x 472 columns
    print('After proportion cutoff: {}'.format(ee.shape))
    original_columns = ee.columns
    ee = ee.apply(func=median_impute,axis=1,result_type='expand')
    ee.columns = original_columns
    cv = ee.apply(func=lambda x:x.values.std()/x.values.mean(),axis=1,result_type='reduce').values
    ee = ee.loc[cv>VARIATION_CUTOFF,:]
    print('After variation cutoff: {}'.format(ee.shape))
    ee.columns = original_columns   # 34683 rows x 472 columns
    ee.T.to_csv(os.path.join(OUTPUT_DIR,'expr.tsv'),sep='\t')


    expr_path = os.path.join(OUTPUT_DIR,'expr.tsv')
    prior_path = K562_PRIOR_NETWORK_PATH
    expr = pd.read_csv(expr_path,sep='\t',index_col=0)  # 472 rows x 34683 columns
    prior = pd.read_csv(prior_path,sep='\t',index_col=0)  # 73157 rows x 221 columns
    common_events = list(set(expr.columns).intersection(set(prior.index)))
    print('number of valid splicing event in inference: {}'.format(len(common_events)))
    print('inference started')
    expr = expr.loc[:,common_events]   # 472 rows x 17836 columns
    prior = prior.loc[common_events,:]  # 17836 rows x 221 columns
    net = prior.stack().reset_index()
    net.columns = ['target','source','weight']
    net = net.loc[net['weight']!=0,:]

    estimate = dc.run_mdt(mat=expr,net=net)
    estimate.to_csv(os.path.join(root_path,'mdt_estimate.txt'),sep='\t')
    print('inference ended, adding metadata')

    # add metadata
    group = pd.read_csv(meta,index_col=0,sep='\t')
    all_meta = {}
    for c in group.columns:
        all_meta[c] = group[c].to_dict()
    estimate = pd.read_csv(os.path.join(root_path,'mdt_estimate.txt'),sep='\t',index_col=0).T
    mi_sample = tuple(estimate.columns.tolist())
    mi_meta = {}
    for k,v in all_meta.items():
        mi_meta[k] = tuple([v[item] for item in mi_sample])
    mi_array = [mi_sample,*mi_meta.values()]
    mi_array_name = ['sample'] + [*mi_meta.keys()]
    mi = pd.MultiIndex.from_arrays(arrays=mi_array,names=mi_array_name)
    estimate.columns = mi
    estimate.to_csv(os.path.join(root_path,'mdt_estimate_morpheus.txt'),sep='\t')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='RNA-SPRINT command line')
    parser.add_argument('--splicing',type=str,default=None,help='path to the splicing EventAnnotation file')
    parser.add_argument('--prop',type=float,default=0.75,help='Only include splicing that present above that proportion cutoff')
    parser.add_argument('--cv',type=float,default=0.1,help='Only include splicing that preserve covariance coeffieint above this cv cutoff')
    parser.add_argument('--prior',type=str,default=None,help='path to the downloaded k562 prior network file')
    parser.add_argument('--meta',type=str,default=None,help='path to a metadata file to facilitate the visualization, two or more column, first is sample name, other columns are different metadata, with informati')
    parser.add_argument('--outdir',type=str,default=None,help='path to output directory')
    args = parser.parse_args()
    main(args)