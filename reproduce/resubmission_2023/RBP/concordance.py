#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os,sys
import subprocess
from tqdm import tqdm
from scipy.stats import hypergeom


def significance(common,reference,query):
    # original, total target is 103582, limma corrected has 103924, let's take the average as the total target M: 103753
    M = 103753
    N = reference.shape[0]
    n = query.shape[0]
    k = common.shape[0]
    p_val = hypergeom.sf(k,M,n,N)
    return p_val

def concordant_analysis(query_file_path,reference_file_path,same_direction=True,pval=0.05):
    # read file
    query = pd.read_csv(query_file_path, sep='\t', index_col=0)
    reference = pd.read_csv(reference_file_path, sep='\t', index_col=0)

    # apply cutoff
    query = query.loc[query['rawp']<pval,:]
    reference = reference.loc[reference['rawp']<pval,:]

    # metric1: absolute overlap, only consider foreground event
    query['foreground'] = [item.split('|')[0] for item in query.index]
    reference['foreground'] = [item.split('|')[0] for item in reference.index]
    common = reference.join(query.set_index('foreground'), on='foreground', lsuffix='_reference', rsuffix='_query',how='inner')

    if not same_direction:
        common['dPSI_query'] = np.negative(common['dPSI_query'].values)


    # metric2: amongst overlap, the fraction of events whose signs agree with each other
    common_dpsi = pd.DataFrame(
        {'foreground': common['foreground'].values, 
         'dpsi_reference': common['dPSI_reference'].values,
         'dpsi_query': common['dPSI_query'].values})
    same_sign = common_dpsi.apply(func=lambda x: True if (x[1] * x[2] > 0) else False, axis=1)
    frac = np.count_nonzero(same_sign.values) / len(same_sign)
    p_val = significance(common,reference,query)
    result = [round(frac,2),round(common.shape[0],2),round(common.shape[0]/query.shape[0],2),round(common.shape[0]/reference.shape[0],2),p_val]
    return result


'''concordance between shRNA K562 profile and high-low SKCM for all examined RBP'''
pwd = os.getcwd()
os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_K562')
sfs = subprocess.run('for file in *.txt; do echo $(basename $file .txt); done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
os.chdir(pwd)
estimate = pd.read_csv('mdt_estimate.txt',sep='\t',index_col=0)
all_sfs = estimate.columns.tolist()
common_sfs = list(set(sfs).intersection(set(all_sfs))) # 191 out of 221
with open('concordance_combine.txt','w') as f:
    f.write('sf\tsign_agreement\tcommon_event\tcommon/query\tcommon/reference\tp_val\n')
    for sf in tqdm(common_sfs):
        query_file_path = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/survival/DAS_folder/Events-dPSI_0.0_rawp/PSI.low_vs_high.txt'
        reference_file_path = os.path.join('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_K562','{}.txt'.format(sf))
        result = concordant_analysis(query_file_path,reference_file_path,False)
        f.write('{}'.format(sf))
        for item in result:
            f.write('\t{}'.format(item))
        f.write('\n')
sys.exit('stop')

# concordance between two ways of doing DAS
# os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output')
# for cat in ['shRNA_K562','shRNA_HepG2','CRISPR_K562','CRISPR_HepG2']:
#     print('processing {}'.format(cat))
#     with open('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/concordance/{}.txt'.format(cat),'w') as f:
#         f.write('sf\tsign_agreement\tcommon_event\tcommon/query\tcommon/reference\tp_val\n')
#         os.chdir(os.path.join(os.getcwd(),'das_within',cat))
#         sfs = subprocess.run('for file in *.txt; do echo $(basename $file .txt); done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#         os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output')
#         for sf in tqdm(sfs):
#             for compare in ['das_within','das_combine']:
#                 os.chdir(os.path.join(os.getcwd(),compare,cat))
#                 if compare == 'das_within':
#                     query_file_path = os.path.join(os.getcwd(),'{}.txt'.format(sf))
#                 elif compare == 'das_combine':
#                     reference_file_path = os.path.join(os.getcwd(),'{}.txt'.format(sf))
#                 os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output')
#             result = concordant_analysis(query_file_path,reference_file_path)
#             f.write('{}\t{}\n'.format(sf,result))
#     os.chdir('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output')


# concordance between ENCODE and AML leucegene cohort
# 1. U2AF1
reference = {'U2AF1_AML_strigent':'/data/salomonis2/NCI-R01/ONCObrowser/Completed/AML-Leucegene/DE_splicing_events/Events-dPSI_0.1_adjp/PSI.Leucegene.U2AF1-Stringent(R1-C2)_vs_Others.txt',
             'U2AF1_AML_extend':'/data/salomonis2/NCI-R01/ONCObrowser/Completed/AML-Leucegene/DE_splicing_events/Events-dPSI_0.1_adjp/PSI.Leucegene.U2AF1-Expanded(R1-C2-extended)_vs_Others.txt',
             'U2AF1_AML_Q157':'/data/salomonis2/NCI-R01/ONCObrowser/Completed/AML-Leucegene/DE_splicing_events/Events-dPSI_0.1_adjp/PSI.Leucegene.U2AF1-Q157_variants_vs_Others.txt',
             'U2AF1_AML_S34':'/data/salomonis2/NCI-R01/ONCObrowser/Completed/AML-Leucegene/DE_splicing_events/Events-dPSI_0.1_adjp/PSI.Leucegene.U2AF1-S34_variants_vs_Others.txt'}

query = {'U2AF1_shRNA_K562_within':'/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/shRNA_K562/U2AF1.txt',
         'U2AF1_shRNA_HepG2_within':'/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/shRNA_HepG2/U2AF1.txt',
         'U2AF1_shRNA_K562_combine':'/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_K562/U2AF1.txt',
         'U2AF1_shRNA_HepG2_combine':'/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_HepG2/U2AF1.txt'}

index_col = []
column_col = []
value_col = []
for qk,qv in query.items():
    for rk,rv in reference.items():
        index_col.append(qk)
        column_col.append(rk)
        value_col.append(concordant_analysis(qv,rv))
final = pd.crosstab(index=np.array(index_col),columns=np.array(column_col),values=np.array(value_col),rownames=['query_file'],colnames=['reference_file'],aggfunc=np.sum)
final.to_csv('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/concordance/U2AF1_concordance.txt',sep='\t')

# 2. HNRNPK
# reference = {'HNRNPK_kd_subtype':'/data/salomonis2/NCI-R01/ONCObrowser/Completed/AML-Leucegene/DE_splicing_events/Events-dPSI_0.1_adjp/PSI.Leucegene.HNRNPK-KD_subtype_vs_Others.txt',
#              'HNRNKP_kd_variant':'/data/salomonis2/NCI-R01/ONCObrowser/Completed/AML-Leucegene/DE_splicing_events/Events-dPSI_0.1_adjp/PSI.Leucegene.HNRNPK-KD_variants_vs_Others.txt'}
# query = {'HNRNPK_shRNA_K562_within':'/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/shRNA_K562/HNRNPK.txt',
#          'HNRNPK_shRNA_HepG2_within':'/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within/shRNA_HepG2/HNRNPK.txt',
#          'HNRNPK_shRNA_K562_combine':'/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_K562/HNRNPK.txt',
#          'HNRNPK_shRNA_HepG2_combine':'/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_combine/shRNA_HepG2/HNRNPK.txt'}    
# index_col = []
# column_col = []
# value_col = []
# for qk,qv in query.items():
#     for rk,rv in reference.items():
#         index_col.append(qk)
#         column_col.append(rk)
#         value_col.append(concordant_analysis(qv,rv))
# final = pd.crosstab(index=np.array(index_col),columns=np.array(column_col),values=np.array(value_col),rownames=['query_file'],colnames=['reference_file'],aggfunc=np.sum)
# final.to_csv('/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/concordance/HNRNPK_concordance.txt',sep='\t')
