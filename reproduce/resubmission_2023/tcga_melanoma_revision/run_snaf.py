#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
sys.path.insert(0,'/data/salomonis2/software')
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import anndata as ad
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import subprocess

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# # # get reduced junction
# df = snaf.get_reduced_junction_matrix(pc='counts.TCGA-SKCM.txt',pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')
# # df.to_csv('../ts/logic_gate/tcga_melanoma_df.txt',sep='\t')

# # run SNAF
# netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
# db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
# tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
# gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
# add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

# snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
# surface.initialize(db_dir=db_dir)


# collect additional neoantigens from the analysis
from Bio.SeqIO.FastaIO import SimpleFastaParser
dic = {}
with open('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/MS_validation_concat_unique2.fasta','r') as handle:
    for title,seq in SimpleFastaParser(handle):
        dic.setdefault(title,[]).append(seq)

os.chdir('/data/salomonis-archive/MS/melanoma/normal_run')
all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
all_patients.remove('Mel-16')
os.chdir('/data/salomonis-archive/MS/melanoma')
data = []
for p in all_patients:
    df = pd.read_csv(os.path.join('/data/salomonis-archive/MS/melanoma/normal_run',p,'combined','txt','peptides.txt'),sep='\t')
    df = df.loc[df['Proteins'].notna(),:]
    df = df.loc[df['Proteins'].str.startswith('>ENSG'),:]
    col = []
    for item1,item2 in zip(df['Sequence'],df['Proteins']):
        if ';' in item2:
            item2 = item2.split(';')[0]
        possible = dic[item2]
        if item1 in possible:
            col.append(True)
        else:
            col.append(False)
    df['exact_same'] = col
    data.append(df)
final = pd.concat(data,axis=0,keys=all_patients).reset_index(level=0)
final.to_csv('additional_neoantigens_strigent.txt',sep='\t',index=None)
sys.exit('stop')

            




# fdr 0.05
db_concat = ['/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot.fasta',
             '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/MS_validation_concat_unique2.fasta']
os.chdir('/data/salomonis-archive/MS/melanoma/normal_run')
all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
for p in all_patients:
    os.chdir('/data/salomonis-archive/MS/melanoma/normal_run/{}'.format(p))
    all_samples = subprocess.run('for file in *.raw; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    inputs = []
    pwd = os.getcwd()
    for s in all_samples:
        inputs.append(os.path.join(pwd,s))
    outdir='/data/salomonis-archive/MS/melanoma/normal_run/{}'.format(p)
    base = '/data/salomonis-archive/MS/ovarian/mqpar.xml'
    snaf.proteomics.set_maxquant_configuration(base=base,dbs=db_concat,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=4,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=11)

# fdr 1
db_concat = ['/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot.fasta',
             '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/MS_validation_concat_unique2.fasta']
os.chdir('/data/salomonis-archive/MS/melanoma/normal_run_FDR1')
all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
for p in all_patients:
    os.chdir('/data/salomonis-archive/MS/melanoma/normal_run_FDR1/{}'.format(p))
    all_samples = subprocess.run('for file in *.raw; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    inputs = []
    pwd = os.getcwd()
    for s in all_samples:
        inputs.append(os.path.join(pwd,s))
    outdir='/data/salomonis-archive/MS/melanoma/normal_run_FDR1/{}'.format(p)
    base = '/data/salomonis-archive/MS/ovarian/mqpar.xml'
    snaf.proteomics.set_maxquant_configuration(base=base,dbs=db_concat,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=4,protein_fdr=1,peptide_fdr=1,site_fdr=1,outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=11)
