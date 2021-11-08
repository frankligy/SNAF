#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
sys.path.append('../')
import snaf
import pandas as pd
import numpy as np
import pickle
from itertools import product

# run SNAF
exon_table = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_Ensembl_exon_add_col.txt'
fasta = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_gene-seq-2000_flank.fa'
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
gtex_db = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/GTEx_junction_counts.h5ad'

snaf.initialize(exon_table=exon_table,fasta=fasta,software_path=netMHCpan_path,gtex_db=gtex_db)

# df = pd.read_csv('./counts.test.txt',sep='\t',index_col=0)
# df.index = [item.split('=')[0] for item in df.index]
# df = df.sample(frac=0.01,axis=0).sample(n=5,axis=1)  # downasample for testing speed purpose
# df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]

# jcmq_schuster = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=5)
# jcmq_schuster.parallelize_run(kind=1)
# sample_to_hla = pd.read_csv('./sample_hla.txt',sep='\t',index_col=2)['HLA'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq_schuster.parallelize_run(kind=4,hlas=hlas)
# jcmq_schuster.serialize(outdir='.',name='after_prediction.p')

jcmq_schuster = snaf.JunctionCountMatrixQuery.deserialize(outdir='.',name='after_prediction.p')
jcmq_schuster.show_neoantigen_burden(outdir='.',name='burden.txt',stage=3,only_peptide=False)
jcmq_schuster.show_neoantigen_frequency(outdir='.',name='frequency.txt',stage=3,only_peptide=True,plot=True,plot_name='frequency.pdf')
jcmq_schuster.show_neoantigen_as_fasta(outdir='.',name='neoantigen.fasta',stage=3)
















