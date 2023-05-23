import os
import sys
import pandas as pd
import numpy as np


# preprocess the dataframe
df = pd.read_csv('/mnt/altanalyze_output/ExpressionInput/counts.original.txt',sep='\t',index_col=0)
df.index = [item.split('=')[0] for item in df.index]
df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]
df.to_csv('/mnt/altanalyze_output/ExpressionInput/counts.original.full.txt',sep='\t')

# filter to EventAnnotation file
ee = [':'.join(item.split('|')[0].split(':')[1:]) for item in pd.read_csv('/mnt/altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t')['UID']]
df = df.loc[df.index.isin(set(ee)),:]
df.to_csv('/mnt/altanalyze_output/ExpressionInput/counts.original.pruned.txt',sep='\t')




