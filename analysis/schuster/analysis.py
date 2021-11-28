#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
import snaf
import pandas as pd
import numpy as np
import pickle


# preprocess the dataframe
df = pd.read_csv('./counts.original.txt',sep='\t',index_col=0)
df.index = [item.split('=')[0] for item in df.index]
df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]

# run SNAF
exon_table = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_Ensembl_exon_add_col.txt'
fasta = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_gene-seq-2000_flank.fa'
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
gtex_db = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/GTEx_junction_counts.h5ad'

# schuster dataset
# paper: https://www.sciencedirect.com/science/article/pii/S1074761317300420?via%3Dihub
snaf.initialize(exon_table=exon_table,fasta=fasta,software_path=netMHCpan_path,gtex_db=gtex_db)
# jcmq_schuster = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=df.shape[1])
# jcmq_schuster.parallelize_run(kind=1)
# sample_to_hla = pd.read_csv('./sample_hla.txt',sep='\t',index_col=2)['HLA'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq_schuster.parallelize_run(kind=4,hlas=hlas)
# jcmq_schuster.serialize(outdir='.',name='after_prediction.p')
# jcmq_schuster = snaf.JunctionCountMatrixQuery.deserialize(outdir='.',name='after_prediction.p')
# jcmq_schuster.show_neoantigen_burden(outdir='.',name='burden_stage3.txt',stage=3,verbosity=1)
# jcmq_schuster.show_neoantigen_frequency(outdir='.',name='frequency_stage3.txt',stage=3,verbosity=1,plot=True,plot_name='frequency_stage3.pdf')
# for sample in jcmq_schuster.junction_count_matrix.columns:
#     print('polishing {} results'.format(sample))
#     jcmq_schuster.show_neoantigen_as_fasta(outdir='./fasta',name='neoantigen_{}.fasta'.format(sample),stage=2,sample=sample)
#     snaf.proteomics.remove_redundant('./fasta/neoantigen_{}.fasta'.format(sample),'./fasta/neoantigen_{}_unique.fasta'.format(sample))
#     snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/chopped_uniprot.fasta',
#                                       fa2_path='./fasta/neoantigen_{}_unique.fasta'.format(sample),outdir='./fasta',
#                                       write_unique2=True,prefix='{}_'.format(sample))


# configure maxquant
## OvCa48
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933726.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#4.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#5.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48/OvCa48_classI_Rep#6.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa48'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

## OvCa53
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933729.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53/OvCa53_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53/OvCa53_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53/OvCa53_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53/OvCa53_classI_Rep#4.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53/OvCa53_classI_Rep#5.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa53'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

## OvCa58
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933728.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58/OvCa58_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58/OvCa58_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58/OvCa58_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58/OvCa58_classI_Rep#4.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58/OvCa58_classI_Rep#5.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa58'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

## OvCa64
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933735.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64/OvCa64_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64/OvCa64_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64/OvCa64_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64/OvCa64_classI_Rep#4.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64/OvCa64_classI_Rep#5.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa64'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

## OvCa65
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933734.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65/OvCa65_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65/OvCa65_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65/OvCa65_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65/OvCa65_classI_Rep#4.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65/OvCa65_classI_Rep#5.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa65'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

## OvCa70
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933738.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70/OvCa70_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70/OvCa70_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70/OvCa70_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70/OvCa70_classI_Rep#4.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70/OvCa70_classI_Rep#5.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa70'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# OvCa80
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5947644.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa80/OvCa80_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa80/OvCa80_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa80/OvCa80_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa80/OvCa80_classI_Rep#4.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa80'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# OvCa84
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5947645.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#4.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#5.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# OvCa99
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5947646.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99/OvCa99_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99/OvCa99_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99/OvCa99_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99/OvCa99_classI_Rep#4.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99/OvCa99_classI_Rep#5.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa99'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# OvCa104
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5947647.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa104/OvCa104_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa104/OvCa104_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa104/OvCa104_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa104/OvCa104_classI_Rep#4.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa104'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# OvCa105
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933743.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105/OvCa105_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105/OvCa105_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105/OvCa105_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105/OvCa105_classI_Rep#4.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105/OvCa105_classI_Rep#5.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa105'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# OvCa109
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933745.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa109/OvCa109_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa109/OvCa109_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa109/OvCa109_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa109/OvCa109_classI_Rep#4.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa109'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# OvCa111
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933736.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa111/OvCa111_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa111/OvCa111_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa111/OvCa111_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa111/OvCa111_classI_Rep#4.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa111'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

# OvCa114
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/fasta/SRR5933737.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114/OvCa114_classI_Rep#1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114/OvCa114_classI_Rep#2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114/OvCa114_classI_Rep#3.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114/OvCa114_classI_Rep#4.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

sys.exit('stop')

# post-inspection and visualization
# row_index = jcmq_schuster.subset.index.tolist().index('ENSG00000170421:E27.5-I27.1')
# col_index = jcmq_schuster.subset.columns.tolist().index('SRR5933738.Aligned.sortedByCoord.out.bed')
# nj = jcmq_schuster.results[col_index][row_index]
# nj.visualize('.','prominent.pdf')
# print(jcmq_schuster.subset.loc['ENSG00000170421:E27.5-I27.1','SRR5933738.Aligned.sortedByCoord.out.bed'])  # 40
# print(snaf.gtex.crude_tumor_specificity('ENSG00000170421:E27.5-I27.1',40))   # 0.17032968997955322
print(snaf.gtex.accurate_tumor_specificity('ENSG00000170421:E27.5-I27.1','mle'))   # 1
print(snaf.gtex.accurate_tumor_specificity('ENSG00000170421:E27.5-I27.1','bayesian'))  # 0.94







