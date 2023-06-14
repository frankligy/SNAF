#!/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/snaf_env/bin/python3.7

import os
import sys
import pandas as pd
import numpy as np
import subprocess

all_file_path = subprocess.run("find ./process/ -type f -name \"abundance_gene.tsv\" -exec echo {} \;",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
lis = []
for item in all_file_path:
    _, _, s, f = item.split('/')
    series = pd.read_csv(item,sep='\t',index_col=0).squeeze()
    series.name = s
    lis.append(series)
df = pd.concat(lis,axis=1)
df.to_csv('gene_tpm.txt',sep='\t')




