#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os

# # dummy 
# pair_info = pd.read_csv('../pair_info_fastq.txt',sep='\t',header=None)
# pair_info.columns = ['pair1','pair2']
# sample_name = [row.pair1 + '_' + row.pair2 + '.Aligned.sortedByCoord.out.bed' for row in pair_info.itertuples(index=False)]
# n_group1 = len(sample_name) // 2
# n_group2 = len(sample_name) - n_group1
# task = 'original'
# pd.DataFrame(data={'sample':sample_name,'groupid':[1]*n_group1 + [2]*n_group2, 'label':['group1']*n_group1 + ['group2']*n_group2}).\
#              to_csv('../altanalyze_output/ExpressionInput/groups.{}.txt'.format(task),sep='\t',header=None,index=None)

# with open('../altanalyze_output/ExpressionInput/comps.{}.txt'.format(task),'w') as f:
#     f.write('1\t2\n')


# # DAS analysis within each batch
# pair_info = pd.read_csv('../pair_info_fastq.txt',sep='\t',header=None)
# pair_info.columns = ['pair1','pair2']
# pair_info_dict = pair_info.set_index(keys='pair1').squeeze().to_dict() 
# meta_k = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
# meta_k.columns = [item.replace(' ','_') for item in meta_k.columns]

# for sf,sub1_df in meta_k.groupby(by='Experiment_target'):
#     for assay, sub2_df in sub1_df.groupby(by='Assay'):
#         for cl, sub3_df in sub2_df.groupby(by='Biosample_term_name'):
#             with open('../das_within/groups.{}.txt'.format('_'.join([sf.split('-')[0],assay.split(' ')[0],cl])),'w') as f:
#                 exp = []
#                 control = []
#                 for row in sub3_df.itertuples():
#                     if row.Paired_end == 1:
#                         exp.append(row.Index + '_' + pair_info_dict[row.Index] + '.Aligned.sortedByCoord.out.bed')
#                         for item in row.Controlled_by.split(','):
#                             c_end = item.split('/')[2]  # one control end
#                             try:
#                                 control.append(c_end + '_' + pair_info_dict[c_end] + '.Aligned.sortedByCoord.out.bed')
#                             except KeyError:
#                                 print(sf,assay,cl,c_end)  # see below
#                                 continue
#                 exp = list(set(exp))
#                 control = list(set(control))
#                 for sample in exp:
#                     f.write('{}\t1\texp\n'.format(sample))
#                 for sample in control:
#                     f.write('{}\t2\tcontrol\n'.format(sample))
            
# '''
# these are because of the mis-typed information in metadata, and they don't affect the group file building, as
# I have manually checked all the control end listed below are included in group file.
# DDX27-human shRNA RNA-seq HepG2 ENCFF278TEH
# EIF3D-human shRNA RNA-seq HepG2 ENCFF201DSF
# LIN28B-human shRNA RNA-seq K562 ENCFF464MHZ
# LIN28B-human shRNA RNA-seq K562 ENCFF023EZN
# '''


# # group file for seperating exp and control
# pair_info = pd.read_csv('../pair_info_fastq.txt',sep='\t',header=None)
# pair_info.columns = ['pair1','pair2']
# pair_info_dict = pair_info.set_index(keys='pair1').squeeze().to_dict() 
# meta_k = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
# meta_k.columns = [item.replace(' ','_') for item in meta_k.columns]
# meta_k = meta_k.loc[meta_k['Paired_end']==1,:]
# exp = []
# control = []
# for row in meta_k.itertuples():
#     exp.append(row.Index + '_' + pair_info_dict[row.Index] + '.Aligned.sortedByCoord.out.bed')
#     for item in row.Controlled_by.split(','):
#         c_end = item.split('/')[2]
#         try:
#             control.append(c_end + '_' + pair_info_dict[c_end] + '.Aligned.sortedByCoord.out.bed')
#         except:
#             print(c_end)
#             continue
# exp = list(set(exp))
# control = list(set(control))
# with open('batch_correction/groups_exp_control.txt','w') as f:
#     for s in exp:
#         f.write('{}\t{}\t{}\n'.format(s,1,'exp'))
#     for s in control:
#         f.write('{}\t{}\t{}\n'.format(s,2,'control'))

# DAS analysis across batch, or combine all controls
pair_info = pd.read_csv('../pair_info_fastq.txt',sep='\t',header=None)
pair_info.columns = ['pair1','pair2']
pair_info_dict = pair_info.set_index(keys='pair1').squeeze().to_dict() 
meta_k = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
meta_k.columns = [item.replace(' ','_') for item in meta_k.columns]

for assay,sub1_df in meta_k.groupby(by='Assay'):
    for cl,sub2_df in sub1_df.groupby(by='Biosample_term_name'):
        # collect all controls
        control = []
        sub2_df = sub2_df.loc[sub2_df['Paired_end']==1,:]
        for item in sub2_df['Controlled_by']:
            for c_end in item.split(','):
                c_end = c_end.split('/')[2]
                try:
                    control.append(c_end + '_' + pair_info_dict[c_end] + '.Aligned.sortedByCoord.out.bed')
                except KeyError:
                    print(assay,cl,c_end)
                    continue
        control = list(set(control))
        # now for each exp
        for sf,sub3_df in sub2_df.groupby(by='Experiment_target'):
            exp = []
            for item in sub3_df.index:
                exp.append(item + '_' + pair_info_dict[item] + '.Aligned.sortedByCoord.out.bed')
            exp = list(set(exp))
            with open('../das_combine/groups.{}.txt'.format('_'.join([sf.split('-')[0],assay.split(' ')[0],cl])),'w') as f:
                for sample in exp:
                    f.write('{}\t1\texp\n'.format(sample))
                for sample in control:
                    f.write('{}\t2\tcontrol\n'.format(sample))                



    















