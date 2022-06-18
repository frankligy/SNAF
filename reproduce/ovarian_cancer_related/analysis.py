import os
import sys
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import pickle
import subprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from copy import deepcopy

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

'''
First, we need to download the fastq files from original PNAS study, and also the MS files as decribed in manuscript method section,
we run STAR 2.6.1 to get the bam files, the command we used for running STAR is as below:

module load STAR/2.6.1# module load STAR/2.6.1
declare -a index_array
index_array+=(SRR5933736 SRR5933737)
for srr in ${index_array[@]}
do
    STAR --genomeDir /data/salomonis2/Genomes/Star-Index-GRCH38/Grch38-STAR-index \
        --readFilesIn $(pwd)/$srr.fastq \
        --outFileNamePrefix $(pwd)/$srr. \
        --runThreadN 4 \
        --outSAMstrandField intronMotif \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile /data/salomonis2/Genomes/Star-Index-GRCH38/Homo_sapiens.GRCh38.85.gtf \
        --limitBAMsortRAM 97417671648
done

Then we can run the first step of AltAnalyze, which has been detailed in the tutorial, after that, you will have a count matrix,
we provide this count matrix in synapse so that we can start the following analysis
The correponding input files and example output files are in https://www.synapse.org/#!Synapse:syn32057190
'''

df = pd.read_csv('count_matrix_ovarian.txt',sep='\t',index_col=0)

# run SNAF-T to get T antigens, all the results will be in /result folder
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'  # REPLACE WITH YOUR NETMHCPAN PATH
db_dir = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'   # REPLACE WITH YOUR DOWNLOADED DATABASE FOLDER

snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='netMHCpan',software_path=netMHCpan_path)
surface.initialize(db_dir=db_dir)

jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=30)
sample_to_hla = pd.read_csv('./sample_hla.txt',sep='\t',index_col=2)['HLA'].to_dict()   # sample_hla.txt is derived from original paper supplemental
hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
jcmq.run(hlas=hlas,outdir='result')
snaf.JunctionCountMatrixQuery.generate_results(path='result/after_prediction.p',outdir='./result')
for sample in jcmq.junction_count_matrix.columns:
    print('polishing {} results'.format(sample))
    jcmq.show_neoantigen_as_fasta(outdir='./fasta',name='neoantigen_{}.fasta'.format(sample),stage=2,verbosity=1,contain_uid=False,sample=sample)
    snaf.proteomics.remove_redundant('./fasta/neoantigen_{}.fasta'.format(sample),'./fasta/neoantigen_{}_unique.fasta'.format(sample))
    snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/human_proteome_uniprot_9_10_mers_unique.fasta', 
                                    # here, human_proteome_uniprot_9_10_mers_unique.fasta is a pre-processed fasta with 9 and 10 mer chopped from uniprot database, 
                                    # this was done by using snaf.proteomics.chop_normal_pep_db function, the original uniprot human proteome fasta (downloaded at 2020 Jan) was provided in synapse
                                      fa2_path='./fasta/neoantigen_{}_unique.fasta'.format(sample),outdir='./fasta',
                                      write_unique2=True,prefix='{}_'.format(sample))

# we need to configure the manquant mqpar.xml files to run the MS validation
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
inputs = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#1.RAW',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#2.RAW',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#3.RAW',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#4.RAW',
          '/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa84/OvCa84_classI_Rep#5.RAW']
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

'''
the command for running maxquant on Linux is

cd /data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114
module load mono/5.20.1
export PATH=$PATH:​/usr/local/mono/5.20.1/bin
mono /data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/MaxQuant_2.0.3.1/bin/MaxQuantCmd.exe /data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/OvCa114/mqpar.xml
'''

# the stacked bar plots and ratio statistics, may need to chagne the relative path a bit, depending on where your MS raw files are
srr_to_id = pd.read_csv('meta.txt',sep='\t',index_col=0).squeeze().to_dict()
fig,ax = plt.subplots()
for i,(srr,id_) in enumerate(srr_to_id.items()):
    db = '{}.Aligned.sortedByCoord.out.bed_unique2.fasta'.format(srr)
    n = subprocess.run(['wc','-l','./fasta/{}'.format(db)],stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[0].split(' ')[0]
    n = int(n)/2
    print(n)
    pep_path = '../../MS/{}/combined/txt/peptides.txt'.format(id_)
    pep = pd.read_csv(pep_path,sep='\t',index_col=0)
    pep = pep.loc[pep['Proteins'].notna(),:]
    v = pep.shape[0]
    v = int(v)
    ax.bar(x=i,height=n,bottom=0,color='b')
    ax.bar(x=i,height=v,bottom=n-v,color='orange')
    ax.text(x=i,y=n+5,s=v,fontsize=5,ha='center')
ax.legend(handles=[Patch(color=i) for i in ['b','orange']],labels=['predicted','MS supported'],loc='upper left',bbox_to_anchor=(1,1),frameon=False)
ax.set_xticks(np.arange(len(srr_to_id)))
ax.set_xticklabels([item.replace('OvCa','P')for item in srr_to_id.values()],rotation=60,fontsize=8)
ax.set_xlabel('Ovarian Cancer Patients')
ax.set_ylabel('Number of Peptides')
plt.savefig('stack_barplot.pdf',bbox_inches='tight')
plt.close()

ms_supported = [352,304,17,81,224,280,462,33,28,207,424,28,267,536]
total = [7112,3653,5186,2405,3200,3116,5312,2804,3524,2928,5023,4263,2128,3256]
rate = np.array(ms_supported) / np.array(total)
'''
[0.04949381 0.08321927 0.00327806 0.03367983 0.07       0.08985879
 0.08697289 0.0117689  0.00794552 0.07069672 0.08441171 0.00656814
 0.12546992 0.16461916]

 0.06342733852873801
'''


# how to generate specific example plots 
uid = 'ENSG00000189180:E2.1-I2.1'    # ENSG00000189180:E2.1-I2.1 # 'ENSG00000170421:E27.1_52900032-E27.3_52899983'  # ENSG00000169908:E2.4-I2.1
print(snaf.uid_to_coord(uid))
snaf.gtex_visual_combine(uid=uid,norm=True,outdir='result',tumor=df,ylim=None)
snaf.gtex_visual_combine(uid=uid,norm=False,outdir='result',tumor=df,ylim=None)
print(snaf.tumor_specificity(uid=uid,method='mean'))
snaf.JunctionCountMatrixQuery.deserialize(name='result/after_prediction.p').visualize(uid,'SRR5947646.Aligned.sortedByCoord.out.bed','result')


# get neoantigen candiates with specificity score, example results provided in synapse
### All candidate
df = pd.read_csv('result/frequency_stage2_verbosity1_uid.txt',sep='\t',index_col=0)
df = snaf.add_gene_symbol_frequency_table(df,remove_quote=True)
df = snaf.add_coord_frequency_table(df,False)
df.to_csv('result/frequency_stage2_verbosity1_uid_gene_coord.txt',sep='\t')
df = pd.read_csv('result/frequency_stage2_verbosity1_uid_gene_coord.txt',sep='\t',index_col=0)
df = snaf.add_tumor_specificity_frequency_table(df,'mean',True)
df = snaf.add_tumor_specificity_frequency_table(df,'mle',False)
df.to_csv('result/frequency_stage2_verbosity1_uid_gene_coord_mean_mle.txt',sep='\t')
df = pd.read_csv('result/frequency_stage2_verbosity1_uid_gene_coord_mean_mle.txt',sep='\t',index_col=0)
df = snaf.add_tumor_specificity_frequency_table(df,'bayesian',True)
df.to_csv('result/frequency_stage2_verbosity1_uid_gene_coord_mean_mle_bayesian.txt',sep='\t')
### sample specific candiate
jcmq = snaf.JunctionCountMatrixQuery.deserialize('result/after_prediction.p')
df = pd.read_csv('result/frequency_stage2_verbosity1_uid_gene_coord_mean_mle.txt',sep='\t',index_col=0)
snaf.report_candidates(jcmq,df,sample='SRR5933726.Aligned.sortedByCoord.out.bed',outdir='result',criterion=[('netMHCpan_el',0,'<=',2),])














