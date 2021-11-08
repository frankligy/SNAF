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

hlas = ['HLA-A*01:01','HLA-A*02:01','HLA-A*24:02','HLA-A*68:01','HLA-B*08:01','HLA-B*08:02']

# start to run
snaf.initialize(exon_table=exon_table,fasta=fasta,software_path=netMHCpan_path,gtex_db=gtex_db)
nj = snaf.NeoJunction(uid='ENSG00000229859:E1.5_61204253-E7.1_61211419',count=25,check_gtex=False)
nj.detect_type()
nj.retrieve_junction_seq()
nj.in_silico_translation()
nj.binding_prediction(hlas=hlas)
nj.immunogenicity_prediction()
for stage,only_peptide in product([1,2,3],[True,False]):
    nj.derive_candidates(stage=stage,only_peptide=only_peptide)
nj.visualize(outdir='.',name='check_nj_visualization.pdf')
for kind in [1,2,3]:
    nj.gtex_viewer(kind=kind,outdir='.')
for method in ['mle','bayesian']:
    nj.infer_tumor_specificity(method=method)

















