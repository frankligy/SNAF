#!/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/snaf_env/bin/python3.7

import os
import sys
import pandas as pd
import numpy as np
import argparse

def main(args):
    path = args.input
    p_dir = os.path.dirname(path)
    df = pd.read_csv(path,sep='\t',index_col=0)
    col1 = []
    col2 = []
    for item in df.index:
        col1.append(item.split('|')[0].split('.')[0])
        col2.append(item.split('|')[1].split('.')[0])
    df['ENST'] = col1
    df['ENSG'] = col2
    df.groupby(by='ENSG')['tpm'].apply(lambda x:x.sum()).to_csv(os.path.join(p_dir,'abundance_gene.tsv'),sep='\t')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert transcript TPM to gene TPM for kallisto output')
    parser.add_argument('--input',type=str,default=None,help='path to the abundance.tsv file')
    args = parser.parse_args()
    main(args)

