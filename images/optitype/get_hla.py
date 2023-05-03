#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import subprocess
import os,sys
import pandas as pd
import numpy as np


samples = pd.read_csv('samples.txt',sep='\t',header=None)
hla_col = []
for s in samples[0]:
    old_cwd = os.getcwd()
    new_cwd = os.path.join(old_cwd,'process',s)
    os.chdir(new_cwd)
    raw_stdout = subprocess.run(['find','.','-type','f','-name','*.tsv'],stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')
    if len(raw_stdout) == 2:
        path_to_result = raw_stdout[0]
        result = pd.read_csv(path_to_result,sep='\t',index_col=0)
        collect = ['HLA-' + hla for hla in result.iloc[0,0:6].tolist()]
        hla_col.append(','.join(collect))
        os.chdir(old_cwd)
    else:
        print(s,raw_stdout)
        os.chdir(old_cwd)
        continue
samples['hla'] = hla_col
samples.rename(columns={0:'sample'},inplace=True)
samples['sample'] = [item+'.R1Aligned.sortedByCoord.out.bed' for item in samples['sample']]
samples.to_csv('sample_hla.txt',sep='\t',index=None)
























