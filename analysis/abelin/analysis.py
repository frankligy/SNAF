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
# hlas = [['HLA-B*57:01'],['HLA-B*51:01'],['HLA-B*54:01'],['HLA-A*29:02']]
# jcmq_schuster.parallelize_run(kind=4,hlas=hlas)
# jcmq_schuster.serialize(outdir='.',name='after_prediction.p')
jcmq_schuster = snaf.JunctionCountMatrixQuery.deserialize(outdir='.',name='after_prediction.p')
# jcmq_schuster.show_neoantigen_burden(outdir='.',name='burden_stage3.txt',stage=3,verbosity=1)
# jcmq_schuster.show_neoantigen_frequency(outdir='.',name='frequency_stage3.txt',stage=3,verbosity=1,plot=True,plot_name='frequency_stage3.pdf')
# for sample in jcmq_schuster.junction_count_matrix.columns:
#     print('polishing {} results'.format(sample))
#     jcmq_schuster.show_neoantigen_as_fasta(outdir='.',name='neoantigen_{}.fasta'.format(sample),stage=2,sample=sample)
#     snaf.proteomics.remove_redundant('neoantigen_{}.fasta'.format(sample),'neoantigen_{}_unique.fasta'.format(sample))
#     snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/chopped_uniprot.fasta',
#                                       fa2_path='./neoantigen_{}_unique.fasta'.format(sample),
#                                       write_unique2=True,prefix='{}_'.format(sample))


# configure maxquant
## A2902
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/RNA/snaf_analysis/fasta/SRR5163127_1_SRR5163127_2.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/A2902/C20150312_HLAIP_1e8ceq_Sonication_2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/A2902/C20150316_HLAIP_1e8ceq_DNase_2.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/A2902'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

## B5101
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/RNA/snaf_analysis/fasta/SRR5163128_1_SRR5163128_2.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5101/C20150312_HLAIP_1e8ceq_Sonication_1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5101/C20150316_HLAIP_1e8ceq_DNase_1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5101/M20150723_CRH_HLA_B51_Rep1_techrep1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5101/M20150723_CRH_HLA_B51_Rep1_techrep2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5101/M20150723_CRH_HLA_B51_Rep2.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5101'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

## B5401
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/RNA/snaf_analysis/fasta/SRR5163129_1_SRR5163129_2.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5401/M20150723_CRH_HLA_B54_Rep1_techrep1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5401/M20150723_CRH_HLA_B54_Rep2_techrep2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5401/C20150414_JGA_HLA_B54_3uLof6uL_5e7ceqoncolumn_rep1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5401/M20150723_CRH_HLA_B54_Rep3_techrep2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5401/C20150414_JGA_HLA_B54_3uLof6uL_5e7ceqoncolumn_rep2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5401/M20150723_CRH_HLA_B54_Rep3_techrep1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5401/M20150723_CRH_HLA_B54_Rep2_techrep1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5401/M20150626_JGA_HLA_B54_biorep2_5e7ceq_reinject_AcOHmobilephase.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5401/M20150626_JGA_HLA_B54_biorep1_5e7ceq_reinject_AcOHmobilephase.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5401'
snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                           outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

## B5701
dbs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/RNA/snaf_analysis/fasta/SRR5163130_1_SRR5163130_2.Aligned.sortedByCoord.out.bed_unique2.fasta']
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5701/C20150414_JGA_HLA_B57_3uLof6uL_5e7ceqoncolumn_rep2.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5701/C20150414_JGA_HLA_B57_star_3uLof6uL_5e7ceqoncolumn_rep1.raw',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5701/M20150626_JGA_HLA_B57_biorep2_5e7ceq_reinject_AcOHmobilephase.raw']
outdir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/abelin/MS/massive.ucsd.edu/MSV000080527/raw/RAW/B5701'
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







