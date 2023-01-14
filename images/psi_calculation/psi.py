import os
import sys
import pandas as pd
import numpy as np
from copy import deepcopy


def find_uid_in_clique(the_uid,region,strand,uid2coords):
    clique = {the_uid:region}
    for uid,coords in uid2coords.items():
        if uid != the_uid:
            if strand == '+':
                cond = (max(int(region[0]),int(coords[0])) - min(int(region[1]),int(coords[1]))) < 0
                if cond:  # they overlap
                    clique[uid] = coords
            elif strand == '-':
                cond = (max(int(region[1]),int(coords[1])) - min(int(region[0]),int(coords[0]))) < 0
                if cond: # they overlap
                    clique[uid] = coords
    return clique

def calculate_psi_core(clique,uid,count,sample_columns):
    sub_count = count.loc[list(clique.keys()),sample_columns]
    if sub_count.shape[0] > 1:
        mat = sub_count.values
        psi = mat[0,:] / mat.sum(axis=0)
        bg_uid = sub_count.index.tolist()[np.argmax(mat[1:,:].sum(axis=1)) + 1]
    else:  # no background event
        mat = sub_count.values
        psi = mat[0,:] / mat.sum(axis=0)
        bg_uid = 'None'      
    return psi,bg_uid
        


def calculate_psi_per_gene(count):
    return_data = []
    count = count.set_index('uid')
    uid2coords = count.apply(lambda x:[x['start'],x['end']],axis=1,result_type='reduce').to_dict()
    for row in count.iterrows():
        uid = row[0]
        region = uid2coords[uid]
        strand = row[1]['strand']
        clique = find_uid_in_clique(uid,region,strand,uid2coords)
        index_list = deepcopy(row[1].index.tolist())
        for c in ['gene','chrom','start','end','strand']:
            index_list.remove(c)
        sample_columns = index_list
        psi_values,bg_uid = calculate_psi_core(clique,uid,count,sample_columns)
        data = (uid,bg_uid,*psi_values)
        return_data.append(data)
    df = pd.DataFrame.from_records(data=return_data,columns=['uid','bg_uid']+sample_columns)
    df.to_csv('check.txt',sep='\t')

        




count = pd.read_csv('counts.original.txt',sep='\t',index_col=0)
col_uid = []
col_gene = []
col_chrom = []
col_start = []  # logical start and end, not physical on the forward strand
col_end = []
col_strand = []
for item in count.index:
    uid,coords = item.split('=')
    chrom, coords = coords.split(':')
    start, end = coords.split('-')
    strand = '+' if start < end else '-'
    gene = uid.split(':')[0]
    col_uid.append(uid)
    col_gene.append(gene)
    col_chrom.append(chrom)
    col_start.append(start)
    col_end.append(end)
    col_strand.append(strand)
for name,col in zip(['uid','gene','chrom','start','end','strand'],[col_uid,col_gene,col_chrom,col_start,col_end,col_strand]):
    count[name] = col

test_count = count.loc[count['gene']=='ENSG00000000003',:]
calculate_psi_per_gene(test_count)

