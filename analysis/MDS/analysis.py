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
gtex_db = '/data/salomonis2/LabFiles/Frank-Li/refactor/data/GTEx_junction_psi.h5ad'
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'

df_raw = pd.read_csv('Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t',index_col='UID').iloc[:,10:]
def median_impute(x):
    med = np.ma.median(np.ma.masked_invalid(x.values))
    result = np.nan_to_num(x.values,nan=med)
    return result
df = df_raw.apply(func=median_impute,axis=1,result_type='expand')
df.columns = df_raw.columns
df.index = [':'.join(item.split('|')[0].split(':')[1:]) for item in df.index]
df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]

meta = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/AML/MDS/MDS-matrix.splicing-metadata.txt',sep='\t',index_col=0)
tmp = meta.loc['healthy',:]
controls = tmp.loc[tmp==1].index
add_control = df.loc[:,controls]   # (56163, 23)

snaf.initialize(exon_table=exon_table,fasta=fasta,gtex_db=gtex_db,software_path=netMHCpan_path,add_control=add_control,binding_method='netMHCpan')
# surface.initialize(transcript_db=transcript_db,exon_table=exon_table,fasta=fasta,
#                        membrane_db=membrane_db,biotype_db=biotype_db,membrane_fasta_db=membrane_fasta_db)

'''B cell neoantigen'''
# df_tumor = df.loc[:,np.logical_not(df.columns.isin(controls))]
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df_tumor,cores=30,add_control=add_control)
# neojunctions = jcmq.valid
# membrane_uid = surface.filter_to_membrane_protein(neojunctions)   # 313 events
# surface.single_run(membrane_uid,2,True,'/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',serialize=True)
# c_candidates,c_further = surface.process_results('single_run_surface_antigen.p',strigency=3)
# print(c_candidates,c_further)

# sa = surface.individual_check('ENSG00000026751:E1.6-E4.1',indices=[3])
# print(sa.full_length)
# print(sa.orfp)

# debug
# df = df.sample(frac=0.01,axis=0).sample(n=5,axis=1)


'''T cell neoantigen'''
# df_tumor = df.loc[:,np.logical_not(df.columns.isin(controls))]
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df_tumor,cores=30,add_control=add_control)
# jcmq.parallelize_run(kind=1)
# sample_to_hla = pd.read_csv('sample_hla.txt',sep='\t',index_col=0)['hla'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq.parallelize_run(kind=3,hlas=hlas)
# jcmq.serialize(outdir='.',name='after_prediction.p')
# jcmq = snaf.JunctionCountMatrixQuery.deserialize(name='./after_prediction.p')
# burden_stage0 = jcmq.cond_df.astype('int8')
# burden_stage0.loc['burden'] = burden_stage0.sum(axis=0).values
# burden_stage0['mean'] = burden_stage0.mean(axis=1).values
# burden_stage0.to_csv('burden_stage0.txt',sep='\t')
# for stage in [3,2,1]:
#     jcmq.show_neoantigen_burden(outdir='.',name='burden_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False)
#     jcmq.show_neoantigen_frequency(outdir='.',name='frequency_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False,plot=True,plot_name='frequency_stage{}.pdf'.format(stage))
# jcmq.show_neoantigen_frequency(outdir='.',name='frequency_stage3_verbosity1_uid.txt',stage=3,verbosity=1,contain_uid=True,plot=False)



# patient analysis
mutation = pd.read_csv('MDS-matrix.splicing-metadata.txt',sep='\t',index_col=0).T
mutation = mutation.stack().to_frame(name='value').reset_index(level=-1)
mutation = mutation.loc[mutation['value']==1,:]   # level_1  value
burden = pd.read_csv('burden_stage3.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]
burden = burden.loc[np.logical_not(burden.index.isin(controls))]
snaf.downstream.mutation_analysis(mode='plot',burden=burden,mutation=mutation,output='patient_analysis/phenotype_stage3.pdf',
                                  n_sample_cutoff=0,gene_column='level_1',genes_to_plot=['SKI-high','RA','RARS','RAEB','RAEB-1','RAEB-2',
                                  'RAEB-T','U2AF1','SF3B1','SRFS2'])





















