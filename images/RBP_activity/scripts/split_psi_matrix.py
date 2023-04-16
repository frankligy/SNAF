#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import os,sys,subprocess
import numpy as np
import pandas as pd


pair_info = pd.read_csv('../pair_info_fastq.txt',sep='\t',header=None)
pair_info.columns = ['pair1','pair2']
pair_info_dict = pair_info.set_index(keys='pair1').squeeze().to_dict() 

meta_k = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
meta_c = pd.read_csv('../processed_metadata_control.tsv',sep='\t',index_col=0)
meta = pd.concat([meta_k,meta_c])
meta.columns = [item.replace(' ','_') for item in meta.columns]
sample_name_by_cat = {}   # {'shRNA_K562':[column_id in eventannotation file,]}
for cl,sub_df in meta.groupby(by='Biosample_term_name'):
    for assay,sub2_df in sub_df.groupby(by='Assay'):
        assay = assay.split(' ')[0]
        for row in sub2_df.itertuples():
            if row.Paired_end == 1:
                file_id = row.Index
                pair_id = pair_info_dict[file_id]
                sample_name = file_id + '_' + pair_id  + '.Aligned.sortedByCoord.out.bed'
                sample_name_by_cat.setdefault('{}_{}'.format(assay,cl),[]).append(sample_name)

psi_raw = pd.read_csv('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t',index_col='UID').iloc[:,10:]
for key,value in sample_name_by_cat.items():
    sub_psi_raw = psi_raw.loc[:,value]
    sub_psi_raw.to_csv('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/AltResults/AlternativeOutput/{}_Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt'.format(key),sep='\t')



