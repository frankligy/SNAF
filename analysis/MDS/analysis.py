#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import pickle



# run SNAF
transcript_db = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/mRNA-ExonIDs.txt'
exon_table = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_Ensembl_exon_add_col.txt'
fasta = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_gene-seq-2000_flank.fa'
membrane_db = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/human_membrane_proteins_acc2ens.txt'
biotype_db = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/Hs_Ensembl_transcript-biotypes.txt'
membrane_fasta_db = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/uniprot_isoform_enhance.fasta'
gtex_db = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/GTEx_junction_counts.h5ad'

df = pd.read_csv('./counts.MDS.txt',sep='\t',index_col=0)
df.index = [item.split('=')[0] for item in df.index]
df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]
ee = [':'.join(item.split('|')[0].split(':')[1:]) for item in pd.read_csv('Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t')['UID']]
df = df.loc[df.index.isin(set(ee)),:]

meta = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/AML/MDS/MDS-matrix.splicing-metadata.txt',sep='\t',index_col=0)
tmp = meta.loc['healthy',:]
controls = tmp.loc[tmp==1].index
add_control = df.loc[:,controls]   # (56163, 23)

snaf.initialize(exon_table=exon_table,fasta=fasta,gtex_db=gtex_db,software_path=None,add_control=add_control)
surface.initialize(transcript_db=transcript_db,exon_table=exon_table,fasta=fasta,
                       membrane_db=membrane_db,biotype_db=biotype_db,membrane_fasta_db=membrane_fasta_db)

'''B cell neoantigen'''
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df)
# neojunctions = jcmq.valid
# membrane_uid = surface.filter_to_membrane_protein(neojunctions)   # 812 events
# surface.single_run(membrane_uid,2,True,'/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',serialize=True)
# c_candidates,c_further = surface.process_results('single_run_surface_antigen.p',strigency=3)
# print(c_candidates,c_further)

'''T cell neoantigen'''
jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=50)
jcmq.parallelize_run(kind=1)
print(jcmq)
sample_to_hla = pd.read_csv('sample_hla.txt',sep='\t',index_col=0)['hla'].to_dict()
hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
jcmq.parallelize_run(kind=3,hlas=hlas)
jcmq.serialize(outdir='.',name='after_prediction.p')
jcmq = snaf.JunctionCountMatrixQuery.deserialize(name='./after_prediction.p')
burden_stage0 = jcmq.cond_df.astype('int8')
burden_stage0.loc['burden'] = burden_stage0.sum(axis=0).values
burden_stage0['mean'] = burden_stage0.mean(axis=1).values
burden_stage0.to_csv('burden_stage0.txt',sep='\t')
for stage in [3,2,1]:
    jcmq.show_neoantigen_burden(outdir='.',name='burden_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False)
    jcmq.show_neoantigen_frequency(outdir='.',name='frequency_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False,plot=True,plot_name='frequency_stage{}.pdf'.format(stage))
jcmq.show_neoantigen_frequency(outdir='.',name='frequency_stage3_verbosity1_uid.txt',stage=3,verbosity=1,contain_uid=True,plot=False)























