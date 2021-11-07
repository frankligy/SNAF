#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
sys.path.append('../../')
import snaf
import pandas as pd
import numpy as np
import pickle


# preprocess the dataframe
df = pd.read_csv('./counts.test.txt',sep='\t',index_col=0)
df.index = [item.split('=')[0] for item in df.index]

# run SNAF
exon_table = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_Ensembl_exon_add_col.txt'
fasta = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_gene-seq-2000_flank.fa'
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
gtex_db = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/GTEx_junction_counts.h5ad'

# schuster dataset
# paper: https://www.pnas.org/content/114/46/E9942.long#sec-9
# GEO: https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=398141
snaf.initialize(exon_table=exon_table,fasta=fasta,software_path=netMHCpan_path,gtex_db=gtex_db)
# for testing purpose, let's downsample the df a bit
df = df.sample(frac=0.01,axis=0).sample(n=5,axis=1)
jcmq_schuster = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=None)
jcmq_schuster.parallelize_run(kind=1)
sample_to_hla = pd.read_csv('./sample_hla.txt',sep='\t',index_col=2)['HLA'].to_dict()
hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
jcmq_schuster.parallelize_run(kind=4,hlas=hlas)
with open('breakpoint.p','wb') as f:
    pickle.dump(jcmq_schuster,f)
with open('breakpoint.p','rb') as f:
    jcmq_schuster = pickle.load(f)
jcmq_schuster.show_neoantigen_burden('./','check_burden.txt',True)











