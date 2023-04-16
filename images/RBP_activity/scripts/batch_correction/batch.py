#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os

import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from functools import reduce
from copy import deepcopy


# construct batch file for each category
count_file = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/ExpressionInput/counts.original.txt'
count = pd.read_csv(count_file,sep='\t',index_col=0)
pair_info = pd.read_csv('../../pair_info_fastq.txt',sep='\t',header=None)
pair_info.columns = ['pair1','pair2']
pair_info_dict = pair_info.set_index(keys='pair1').squeeze().to_dict() 
meta_k = pd.read_csv('../../processed_metadata_KD.tsv',sep='\t',index_col=0)
meta_k.columns = [item.replace(' ','_') for item in meta_k.columns]
for assay,sub1_df in meta_k.groupby(by='Assay'):
    for cl,sub2_df in sub1_df.groupby(by='Biosample_term_name'):
        # for each category, only look at pair1,becasue pair2 will be auto-filled by pair_info_dict
        sub2_df = sub2_df.loc[sub2_df['Paired_end']==1,:]
        batch_partition = []  # each item is a container, each container can include overlapping keys (key means unique control_by column value, but as list split by ,)
        for cb, sub3_df in sub2_df.groupby(by='Controlled_by'):
            cb = cb.split(',')
            is_cb_overlap = False
            for i,c in enumerate(batch_partition): 
                for k in c:
                    if len(set(cb).intersection(set(k))) > 0:
                        is_cb_overlap = True
                        batch_partition[i].append(cb)
                        break
            if not is_cb_overlap:   # this cb doesn't overlap with any existing batch
                batch_partition.append([cb])
        batch_dict = {}  # key is a unique control_by column value, value is batchid
        for i,c in enumerate(batch_partition):
            for k in c:
                k = ','.join(k)
                batch_dict[k] = 'batch{}'.format(i+1)
        with open('batch_{}_{}.txt'.format(assay.split(' ')[0],cl),'w') as f:
            for row in sub2_df.itertuples():
                sample = row.Index + '_' + pair_info_dict[row.Index] + '.Aligned.sortedByCoord.out.bed'
                cb_key = row.Controlled_by
                batch_id = batch_dict[cb_key]
                f.write('{}\t{}\n'.format(sample,batch_id))
            queried_ends = []
            for k,v in batch_dict.items():
                for kk in k.split(','):
                    c_end = kk.split('/')[2]
                    if c_end not in queried_ends:
                        try:
                            sample = c_end + '_' + pair_info_dict[c_end] + '.Aligned.sortedByCoord.out.bed'
                        except:
                            print(assay,cl,c_end)  # see note below
                            continue
                        batch_id = v
                        queried_ends.append(c_end)
                        f.write('{}\t{}\n'.format(sample,batch_id))
        # make sure the order is the same as the count file, so the r can read
        batch = pd.read_csv('batch_{}_{}.txt'.format(assay.split(' ')[0],cl),sep='\t',index_col=0,header=None)
        batch = batch.squeeze().to_dict()
        count_s = count.loc[:,list(batch.keys())]
        correspond_batch = [item.split('batch')[1] for item in count_s.columns.map(batch).values]
        final = pd.DataFrame({'file':count_s.columns.values,'batch':correspond_batch})
        final.to_csv('batch_{}_{}.txt'.format(assay.split(' ')[0],cl),sep='\t',index=None)



'''
Again, four mis-typed control samples
shRNA RNA-seq HepG2 ENCFF278TEH
shRNA RNA-seq HepG2 ENCFF201DSF
shRNA RNA-seq K562 ENCFF023EZN
shRNA RNA-seq K562 ENCFF464MHZ
'''