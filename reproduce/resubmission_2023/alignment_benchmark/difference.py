#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os,sys
# sys.path.insert(0,'/data/salomonis2/software')
# import snaf
# from snaf import surface
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import subprocess
from matplotlib_venn import venn2

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank')

'''
TCGA-A1-A0SI-01A-11R-A144-07
TCGA-A1-A0SM-01A-11R-A084-07
'''


# construct
'''original bam'''
df1 = pd.read_csv('/data/salomonis2/NCI-R01/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/ExpressionInput/counts.TCGA-BRCA.txt',sep='\t',
                 index_col=0,usecols=['AltAnalyze_ID','TCGA-A1-A0SI-01A-11R-A144-07.bed','TCGA-A1-A0SM-01A-11R-A084-07.bed'])

'''STAR 2.6.1 one pass, ensembl 85'''
df2 = pd.read_csv('/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/TCGA-BRCA-STAR/run/altanalyze_output/ExpressionInput/counts.original.txt',sep='\t',index_col=0)
df2.rename(columns=lambda x:'.'.join([x.split('.')[0],x.split('.')[-1]]),inplace=True)

'''STAR 2.6.1 intermediate using all SJ out, Gencode v36'''
df3 = pd.read_csv('/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/TCGA-BRCA-STAR2pass-GDCRef/bams/ExpressionInput/counts.test.txt',sep='\t',index_col=0)
df3.rename(columns=lambda x:'.'.join([x.split('.')[0],x.split('.')[-1]]),inplace=True)

'''STAR 2.6.1 full two pass, Gencode v36'''
df4 = pd.read_csv('/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank/bam/altanalyze_output/ExpressionInput/counts.original.txt',sep='\t',index_col=0)
df4.rename(columns=lambda x:'.'.join([x.split('_')[0],x.split('.')[-1]]),inplace=True)

'''STAR 2.4.0 full two pass, Gencode v36'''
df5 = pd.read_csv('/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank_second_try/run/altanalyze_output/ExpressionInput/counts.original.txt',sep='\t',index_col=0)
df5.rename(columns=lambda x:'.'.join([x.split('_')[0],x.split('.')[-1]]),inplace=True)

'''STAR 2.4.0 intermediate using all SJ out, Gencode v36'''
df6 = pd.read_csv('/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/STAR2.4.0/bams/ExpressionInput/counts.test.txt',sep='\t',index_col=0)
df6.rename(columns=lambda x:'.'.join([x.split('.')[0],x.split('.')[-1]]),inplace=True)

'''STAR 2.7.5 intermediate using all SJ out, Gencode v36'''
df7 = pd.read_csv('/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/STAR2.7.5c/bams/ExpressionInput/counts.test.txt',sep='\t',index_col=0)
df7.rename(columns=lambda x:'.'.join([x.split('.')[0],x.split('.')[-1]]),inplace=True)


# test
sample = 'TCGA-A1-A0SI-01A-11R-A144-07.bed'
cutoff_e = 10
cutoff_t = 0.1  # percentage

fractions = []
base_series = df1[sample]
base_series = base_series.loc[base_series>cutoff_e]
fig,axes = plt.subplots(ncols=4,nrows=1,gridspec_kw={'wspace':0.5},figsize=(24,4.8))
for i,(df,id_) in enumerate(zip([df2,df3,df4,df5],['2.6.1\none pass\nEnsembl 85','2.6.1\nintermediate two pass\nGencode v36','2.6.1\nfull two pass\nGencode v36','2.4.0\nfull two pass\nGencode v36'])):
    series = df[sample]
    series = series.loc[series>cutoff_e]
    venn2(subsets=[set(base_series.index),set(series.index)],set_labels=['TCGA BAM',id_],ax=axes[i])
    df = pd.concat([base_series,series],axis=1,join='inner',ignore_index=True)
    abs_diff = df.apply(func=lambda x:abs(x.iloc[0] - x.iloc[1]),axis=1,result_type='reduce').values
    tolerable_diff = df.apply(func=lambda x:x.max() * cutoff_t,axis=1,result_type='reduce').values
    cond = abs_diff <= tolerable_diff
    df['abs_diff'] = abs_diff
    df['tolerable_diff'] = tolerable_diff
    df['cond'] = cond
    fraction = round(df['cond'].sum() / df.shape[0],4)
    fractions.append(fraction)
    axes[i].set_title('fraction: {}'.format(fraction))
fig.suptitle('cutoff_e: {} cutoff_t: {}'.format(cutoff_e,cutoff_t))
plt.savefig('/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank_second_try/difference/cutoff_e_{}_cutoff_t_{}.pdf'.format(cutoff_e,cutoff_t),bbox_inches='tight')
plt.close()

fig,ax = plt.subplots()
ax.bar(x=[1,2,3,4],height=fractions)
plt.savefig('/data/salomonis-archive/FASTQs/NCI-R01/TestRun_Aligner/Frank_second_try/difference/fractions.pdf',bbox_inches='tight')
plt.close()








