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
from copy import deepcopy
import anndata as ad

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# # get junction count matrix
# df = snaf.get_reduced_junction_matrix(pc='counts.original.txt',pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')

# # run SNAF
# netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
# db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
# tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
# gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
# add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

# snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
# surface.initialize(db_dir=db_dir)

# # T antigen
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=30,add_control=add_control,outdir='result',filter_mode='maxmin')
# sample_to_hla = pd.read_csv('./sample_hla.txt',sep='\t',index_col=2)['HLA'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq.run(hlas=hlas,outdir='./result')
# snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')
# # proteome
# jcmq = snaf.JunctionCountMatrixQuery.deserialize('./result/after_prediction.p')
# for sample in jcmq.junction_count_matrix.columns:
#     print('polishing {} results'.format(sample))
#     jcmq.show_neoantigen_as_fasta(outdir='./fasta',name='neoantigen_{}.fasta'.format(sample),stage=2,verbosity=1,contain_uid=False,sample=sample)
#     snaf.proteomics.remove_redundant('./fasta/neoantigen_{}.fasta'.format(sample),'./fasta/neoantigen_{}_unique.fasta'.format(sample))
#     snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot_9_10_mers_unique.fasta',
#                                       fa2_path='./fasta/neoantigen_{}_unique.fasta'.format(sample),outdir='./fasta',
#                                       write_unique2=True,prefix='{}_'.format(sample))



# configure maxquant
fasta_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ovarian/fasta'
raw_dir = '/data/salomonis-archive/MS/ovarian'
root_dir = os.path.dirname(os.path.abspath(__file__))
dic = {
    'OvCa48':'SRR5933726',
    'OvCa53':'SRR5933729',
    'OvCa58':'SRR5933728',
    'OvCa64':'SRR5933735',
    'OvCa65':'SRR5933734',
    'OvCa70':'SRR5933738',
    'OvCa80':'SRR5947644',
    'OvCa84':'SRR5947645',
    'OvCa99':'SRR5947646',
    'OvCa104':'SRR5947647',
    'OvCa105':'SRR5933743',
    'OvCa109':'SRR5933745',
    'OvCa111':'SRR5933736',
    'OvCa114':'SRR5933737'
}

# for k,v in dic.items():
#     # both, fdr 0.05
#     dbs = ['/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot.fasta',
#            os.path.join(fasta_dir,'{}_secondAligned.sortedByCoord.out.bed_unique2.fasta'.format(v))]
#     os.chdir(os.path.join(raw_dir,k))
#     inputs = subprocess.run("for file in *.raw; do echo $file; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     os.chdir(root_dir)
#     inputs = [os.path.join(raw_dir,k,inp) for inp in inputs]
#     outdir = os.path.join(raw_dir,k)
#     base = '/data/salomonis-archive/MS/ovarian/mqpar.xml'
#     snaf.proteomics.set_maxquant_configuration(base=base,dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=4,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=11)
#     # both, fdr 1
#     dbs = ['/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot.fasta',
#            os.path.join(fasta_dir,'{}_secondAligned.sortedByCoord.out.bed_unique2.fasta'.format(v))]
#     os.chdir(os.path.join(raw_dir,'{}_FDR1'.format(k)))
#     inputs = subprocess.run("for file in *.raw; do echo $file; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     os.chdir(root_dir)
#     inputs = [os.path.join(raw_dir,'{}_FDR1'.format(k),inp) for inp in inputs]
#     outdir = os.path.join(raw_dir,'{}_FDR1'.format(k))
#     base = '/data/salomonis-archive/MS/ovarian/mqpar.xml'
#     snaf.proteomics.set_maxquant_configuration(base=base,dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=4,protein_fdr=1,peptide_fdr=1,site_fdr=1,outdir=outdir,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=11)


# analyze
def is_neoantigen(s):
    b = (~s.str.startswith('REV_')) & (~s.str.startswith('CON_')) & (~s.str.startswith('sp|')) & (s!=' ')
    return b

# test
db1_data = []  # for peptide
db2_data = []  # for PSM
 

for sample in dic.keys():
    df1 = pd.read_csv('/data/salomonis-archive/MS/ovarian/{}/combined/txt/msmsScans.txt'.format(sample),sep='\t')
    df2 = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ovarian/MS/{}/combined/txt/msmsScans.txt'.format(sample),sep='\t') 
    df3 = pd.read_csv('/data/salomonis-archive/MS/ovarian/{}_FDR1/combined/txt/msmsScans.txt'.format(sample),sep='\t')
    result_list = []
    for df in [df1,df2,df3]:
        df['id'] = [item1 + ',' + str(item2) for item1,item2 in zip(df['Raw file'],df['Scan number'])]
        df.set_index(keys='id')
        df = df.loc[:,['Raw file','Scan number','Identified','Sequence','Score','PEP','Proteins']]
        result_list.append(df)
    final = pd.concat(result_list,axis=1,keys=['both_FDR0.05','neoantigen_FDR0.05','both_FDR1'])
    final.columns = [l1 + ',' + l2 for l1,l2 in final.columns.tolist()]
    final['is_neoantigen_1'] = is_neoantigen(final['both_FDR0.05,Proteins']) 
    final['is_neoantigen_2'] = is_neoantigen(final['neoantigen_FDR0.05,Proteins'])
    final['is_neoantigen_3'] = is_neoantigen(final['both_FDR1,Proteins'])

    final['identified_type1'] = (final['is_neoantigen_1']) & (final['both_FDR0.05,Identified'] == '+') & (final['both_FDR0.05,Sequence'].values == np.array([item.split('|')[0] for item in final['both_FDR0.05,Proteins']]))
    final['identified_type2'] = (final['neoantigen_FDR0.05,Identified'] == '+') & (final['neoantigen_FDR0.05,Score'] >= 40) & (final['both_FDR0.05,Identified'] != '+')
    final['identified_type3'] = final['is_neoantigen_3'] & (final['both_FDR1,Sequence'].values == np.array([item.split('|')[0] for item in final['both_FDR1,Proteins']]))

    final['additional_PSM'] = (final['neoantigen_FDR0.05,Identified'] == '+') & (final['both_FDR0.05,Identified'] != '+')
    final['additional_PSM_FDR1'] = (final['neoantigen_FDR0.05,Identified'] == '+') & (final['both_FDR0.05,Identified'] != '+') & (final['both_FDR1,Identified'] != '+')

    # calculate number of peptides
    n1 = len(set(final.loc[final['identified_type1'],:]['both_FDR0.05,Sequence'].tolist()))
    n2 = len(set(final.loc[final['identified_type2'],:]['neoantigen_FDR0.05,Sequence'].tolist()))
    n3 = len(set(final.loc[final['identified_type3'],:]['both_FDR1,Sequence'].tolist()))
    total = n1 + n2 + n3
    db1_data.append((sample,total))

    # calculate number of additional PSM
    n1 = final.loc[final['additional_PSM'],:].shape[0]
    n2 = final.loc[final['additional_PSM_FDR1'],:].shape[0]
    db2_data.append((sample,n1,n2))
    # final.to_csv('check_{}.txt'.format(sample),sep='\t')

db1 = pd.DataFrame.from_records(db1_data,columns=['sample','number_of_neoantigens'])
db2 = pd.DataFrame.from_records(db2_data,columns=['sample','additional_PSM_FDR0.05','additional_PSM_FDR1'])

db1.to_csv('db1.txt',sep='\t',index=None)
db2.to_csv('db2.txt',sep='\t',index=None)