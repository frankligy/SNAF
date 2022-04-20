#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm


# preprocess the dataframe
df = pd.read_csv('./counts.TCGA-SKCM.txt',sep='\t',index_col=0)
df.index = [item.split('=')[0] for item in df.index]
df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]

# filter to EventAnnotation file
ee = [':'.join(item.split('|')[0].split(':')[1:]) for item in pd.read_csv('Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t')['UID']]
df = df.loc[df.index.isin(set(ee)),:]


# run SNAF
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
db_dir = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'

snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='netMHCpan',software_path=netMHCpan_path)
surface.initialize(db_dir=db_dir)


'''B cell neoantigen'''
# membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df)
# surface.run(membrane_tuples,outdir='result',tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm')
# surface.generate_results(pickle_path='./result/surface_antigen.p',outdir='result',strigency=5,gtf='./SQANTI-all/collapse_isoforms_classification.filtered_lite.gtf')   # './SQANTI-all/collapse_isoforms_classification.filtered_lite.gtf'
surface.run_dash_B_antigen(pkl='result/surface_antigen.p',candidates='result/candidates_5.txt',python_executable='/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7')
sys.exit('stop')


'''T cell neoantigen'''
#### Step1: Running the program
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=50)
# sample_to_hla = pd.read_csv('../HLA_genotype/sample_hla.txt',sep='\t',index_col=0)['hla'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq.run(hlas=hlas,outdir='./result')
# snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')


### Step2: necessary secondary results
# df = pd.read_csv('./result/frequency_stage3_verbosity1_uid.txt',sep='\t',index_col=0)
# snaf.downstream.add_gene_symbol_frequency_table(df=df).to_csv('./result/frequency_stage3_verbosity1_uid_gene_symbol.txt',sep='\t')
# jcmq = snaf.JunctionCountMatrixQuery.deserialize(name='./result/after_prediction.p')
# jcmq.visualize(uid='ENSG00000167291:E38.6-E39.1',sample='TCGA-DA-A1I1-06A-12R-A18U-07.bed',outdir='./result')

### Step3: downstream analysis (patient level and neoantigen level)

'''patient analysis'''
# 1. survival analysis
# survival = pd.read_csv('patient_analysis/TCGA-SKCM.survival.tsv',sep='\t',index_col=0)  # 463
# burden = pd.read_csv('burden_stage2.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
# burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
# burden,burden_encode,burden_vc = snaf.survival_analysis(burden,survival,n=2,stratification_plot='patient_analysis/stage2_stratify.pdf',survival_plot='patient_analysis/stage2_survival.pdf')
# burden_encode.to_csv('to_anu/stage2_stratify.txt',sep='\t',header=None)

# 2. mutation analysis
# mutation = pd.read_csv('patient_analysis/TCGA-SKCM.mutect2_snv.tsv',sep='\t',index_col=0)  # 467 samples have mutations
# mutation = mutation.loc[mutation['filter']=='PASS',:]
# burden = pd.read_csv('burden_stage3.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
# burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
# # snaf.mutation_analysis(mode='compute',burden=burden,mutation=mutation,output='patient_analysis/stage3_mutation.txt')
# snaf.mutation_analysis(mode='plot',burden=burden,mutation=mutation,output='patient_analysis/stage3_mutation_CAMKK2.pdf',genes_to_plot=['CAMKK2'])


'''neoantigen analysis'''
# 1. physicalchemical properties and occurence frequency
# snaf.downstream.analyze_neoantigens(freq_path='result/frequency_stage2_verbosity1_uid.txt',junction_path='result/burden_stage0.txt',total_samples=472,outdir='result',mers=[9,10],fasta=True)
# snaf.run_dash_T_antigen(input_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/result/shared_vs_unique_neoantigen_all.txt')

        




























