#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
sys.path.insert(0,'/data/salomonis2/software')
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import pickle
import subprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# preprocess the dataframe
df = pd.read_csv('./counts.original.txt',sep='\t',index_col=0)
df.index = [item.split('=')[0] for item in df.index]
df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]

# run SNAF
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
db_dir = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'

snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='netMHCpan',software_path=netMHCpan_path)
surface.initialize(db_dir=db_dir)


# schuster dataset
# paper: 
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=30)
# sample_to_hla = pd.read_csv('./sample_hla.txt',sep='\t',index_col=2)['HLA'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq.run(hlas=hlas,outdir='result')
# snaf.JunctionCountMatrixQuery.generate_results(path='result/after_prediction.p',outdir='./result')
# for sample in jcmq.junction_count_matrix.columns:
#     print('polishing {} results'.format(sample))
#     jcmq.show_neoantigen_as_fasta(outdir='./fasta',name='neoantigen_{}.fasta'.format(sample),stage=2,verbosity=1,contain_uid=False,sample=sample)
#     snaf.proteomics.remove_redundant('./fasta/neoantigen_{}.fasta'.format(sample),'./fasta/neoantigen_{}_unique.fasta'.format(sample))
#     snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/chopped_uniprot.fasta',
#                                       fa2_path='./fasta/neoantigen_{}_unique.fasta'.format(sample),outdir='./fasta',
#                                       write_unique2=True,prefix='{}_'.format(sample))


# df = pd.read_csv('result/frequency_stage2_verbosity1_uid.txt',sep='\t',index_col=0)
# df = snaf.add_gene_symbol_frequency_table(df,remove_quote=True)
# df = snaf.add_coord_frequency_table(df,False)
# df.to_csv('result/frequency_stage2_verbosity1_uid_gene_coord.txt',sep='\t')
# df = pd.read_csv('result/frequency_stage2_verbosity1_uid_gene_coord.txt',sep='\t',index_col=0)
# df = snaf.add_tumor_specificity_frequency_table(df,'mean',True)
# df = snaf.add_tumor_specificity_frequency_table(df,'mle',False)
# df.to_csv('result/frequency_stage2_verbosity1_uid_gene_coord_mean_mle.txt',sep='\t')
df = pd.read_csv('result/frequency_stage2_verbosity1_uid_gene_coord_mean_mle.txt',sep='\t',index_col=0)
df = snaf.add_tumor_specificity_frequency_table(df,'bayesian',True)
df.to_csv('result/frequency_stage2_verbosity1_uid_gene_coord_mean_mle_bayesian.txt',sep='\t')
sys.exit('stop')


# # configure maxquant
# ## OvCa48
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933726.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#4.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#5.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#6.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# ## OvCa53
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933729.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53/OvCa53_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53/OvCa53_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53/OvCa53_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53/OvCa53_classI_Rep#4.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53/OvCa53_classI_Rep#5.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# ## OvCa58
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933728.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58/OvCa58_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58/OvCa58_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58/OvCa58_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58/OvCa58_classI_Rep#4.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58/OvCa58_classI_Rep#5.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# ## OvCa64
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933735.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64/OvCa64_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64/OvCa64_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64/OvCa64_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64/OvCa64_classI_Rep#4.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64/OvCa64_classI_Rep#5.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# ## OvCa65
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933734.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65/OvCa65_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65/OvCa65_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65/OvCa65_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65/OvCa65_classI_Rep#4.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65/OvCa65_classI_Rep#5.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# ## OvCa70
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933738.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70/OvCa70_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70/OvCa70_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70/OvCa70_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70/OvCa70_classI_Rep#4.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70/OvCa70_classI_Rep#5.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# # OvCa80
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5947644.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa80/OvCa80_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa80/OvCa80_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa80/OvCa80_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa80/OvCa80_classI_Rep#4.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa80'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# # OvCa84
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5947645.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#1.RAW',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#2.RAW',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#3.RAW',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#4.RAW',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#5.RAW']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# # OvCa99
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5947646.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99/OvCa99_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99/OvCa99_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99/OvCa99_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99/OvCa99_classI_Rep#4.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99/OvCa99_classI_Rep#5.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# # OvCa104
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5947647.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa104/OvCa104_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa104/OvCa104_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa104/OvCa104_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa104/OvCa104_classI_Rep#4.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa104'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# # OvCa105
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933743.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105/OvCa105_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105/OvCa105_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105/OvCa105_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105/OvCa105_classI_Rep#4.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105/OvCa105_classI_Rep#5.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# # OvCa109
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933745.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa109/OvCa109_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa109/OvCa109_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa109/OvCa109_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa109/OvCa109_classI_Rep#4.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa109'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# # OvCa111
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933736.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa111/OvCa111_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa111/OvCa111_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa111/OvCa111_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa111/OvCa111_classI_Rep#4.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa111'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# # OvCa114
# dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933737.Aligned.sortedByCoord.out.bed_unique2.fasta']
# inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114/OvCa114_classI_Rep#1.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114/OvCa114_classI_Rep#2.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114/OvCa114_classI_Rep#3.raw',
#           '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114/OvCa114_classI_Rep#4.raw']
# outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114'
# snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                            outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)


# stack barplot
# srr_to_id = pd.read_csv('../meta.txt',sep='\t',index_col=0).squeeze().to_dict()
# fig,ax = plt.subplots()
# for i,(srr,id_) in enumerate(srr_to_id.items()):
#     db = '{}.Aligned.sortedByCoord.out.bed_unique2.fasta'.format(srr)
#     n = subprocess.run(['wc','-l','./fasta/{}'.format(db)],stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[0].split(' ')[0]
#     n = int(n)/2
#     print(n)
#     pep_path = '../../MS/{}/combined/txt/peptides.txt'.format(id_)
#     pep = pd.read_csv(pep_path,sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     v = pep.shape[0]
#     v = int(v)
#     ax.bar(x=i,height=n,bottom=0,color='b')
#     ax.bar(x=i,height=v,bottom=n-v,color='orange')
#     ax.text(x=i,y=n+5,s=v,fontsize=5,ha='center')
# ax.legend(handles=[Patch(color=i) for i in ['b','orange']],labels=['predicted','MS supported'],loc='upper left',bbox_to_anchor=(1,1),frameon=False)
# ax.set_xticks(np.arange(len(srr_to_id)))
# ax.set_xticklabels([item.replace('OvCa','P')for item in srr_to_id.values()],rotation=60,fontsize=8)
# ax.set_xlabel('Ovarian Cancer Patients')
# ax.set_ylabel('Number of Peptides')
# plt.savefig('stack_barplot.pdf',bbox_inches='tight')
# plt.close()
    
ms_supported = [352,304,17,81,224,280,462,33,28,207,424,28,267,536]
total = [7112,3653,5186,2405,3200,3116,5312,2804,3524,2928,5023,4263,2128,3256]
rate = np.array(ms_supported) / np.array(total)
'''
[0.04949381 0.08321927 0.00327806 0.03367983 0.07       0.08985879
 0.08697289 0.0117689  0.00794552 0.07069672 0.08441171 0.00656814
 0.12546992 0.16461916]

 0.06342733852873801
'''

# /Volumes/salomonis-archive/BAMs/PublicDatasets/E-MTAB-2836-Grch38_Deep-Healthy-PanTissue/ERR315460_1.bam


# FOR figure2 specific examples
uid = 'ENSG00000170421:E27.1_52900032-E27.3_52899983'    # ENSG00000189180:E2.1-I2.1
print(snaf.uid_to_coord(uid))
snaf.gtex_visual_combine(uid=uid,norm=True,outdir='result',tumor=df,ylim=None)
snaf.gtex_visual_combine(uid=uid,norm=False,outdir='result',tumor=df,ylim=None)
print(snaf.tumor_specificity(uid=uid,method='mean'))
snaf.JunctionCountMatrixQuery.deserialize(name='result/after_prediction.p').visualize(uid,'SRR5947645.Aligned.sortedByCoord.out.bed','result')











