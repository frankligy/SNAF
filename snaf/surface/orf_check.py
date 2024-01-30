import pandas as pd
import numpy as np
import os,sys
import subprocess
import multiprocessing
import re
# for biopython, pip install biopython
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq


def set_global_env(df_exonlist_arg,dict_exonCoords_arg,dict_fa_arg,dict_biotype_arg):
    global df_exonlist
    global dict_exonCoords
    global dict_fa
    global dict_biotype
    df_exonlist = df_exonlist_arg
    dict_exonCoords = dict_exonCoords_arg
    dict_fa = dict_fa_arg
    dict_biotype = dict_biotype_arg

def query_from_dict_fa(abs_start,abs_end,EnsID,strand):
    '''
    abs_start and abs_end always means the xth base in forward strand
    the returned exon_seq, however, means the 5'-3' seq depending on the strand information.
    '''
    if strand == '+':        
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq = dict_fa[EnsID][3]
        start_index = int(abs_start) - start + 2000
        end_index = int(abs_end) - start + 1 + 2000
        exon_seq = seq[start_index:end_index]
    
    elif strand == '-':
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq_reverse = dict_fa[EnsID][3]
        seq_forward = str(Seq(seq_reverse).reverse_complement())  # Hs_gene.fa restore the reverse strand info
        start_index = int(abs_start) - start + 2000
        end_index = int(abs_end) - start + 1 + 2000 # endpoint in python is non-inclusive
        exon_seq_1 = seq_forward[start_index:end_index]
        s = Seq(exon_seq_1)
        exon_seq = str(s.reverse_complement())
    return exon_seq


def get_exon_sequence(exon,ensgid):
    attrs = dict_exonCoords[ensgid][exon]
    return query_from_dict_fa(attrs[2],attrs[3],ensgid,attrs[1])


def back_traverse_exonlist(exonlist,ensgid,n_stride):
    exonlist = exonlist.split('|')
    exon_iterator = iter(exonlist[::-1])
    current_subexon = next(exon_iterator,'end')  # arrive the last one
    current_exon,current_index = current_subexon[1:].split('.')   # str(15)  # str(1)
    current_exon,current_index = int(current_exon),int(current_index) # int(15)  # int(1)
    current_subexon_seq = get_exon_sequence(current_subexon,ensgid)
    n_traversed_exon = 0
    traversed_bases = len(current_subexon_seq)
    while True:
        if n_traversed_exon >= n_stride:
            break
        subexon = next(exon_iterator,'end')
        current_subexon = subexon
        if subexon == 'end':
            break
        else:
            exon,index = subexon[1:].split('.')
            exon,index = int(exon),int(index)
            subexon_seq = get_exon_sequence(subexon,ensgid)
            traversed_bases += len(subexon_seq)
            if exon < current_exon or index < current_index - 1:   # from E15.1 to E14.2 or from E14.5 to E14.3
                n_traversed_exon += 1
                current_exon = exon
                current_index = index
    return traversed_bases,n_traversed_exon

def nmd_check(uid,full_length,orf,n_stride):
    check = []
    ensgid = uid.split(':')[0]
    df_certain = df_exonlist.loc[df_exonlist['EnsGID']==ensgid,['EnsPID','Exons']]
    for o,t,el in zip(orf,full_length,df_certain['Exons']):
        if o != 'unrecoverable' and o != '':
            orf_last_base_pos = t.index(o) + (len(o) - 1)    # think about the start pos is 2, len is 3, means forward 2 (3-1) position
            relative_from_end = len(t) - orf_last_base_pos   # draw a simple graph and think about that.
            traversed_bases, n_traversed_exon = back_traverse_exonlist(el,ensgid,n_stride)
            if n_traversed_exon == n_stride:
                if traversed_bases < relative_from_end:
                    check.append('*')
                else:
                    check.append('#')
            else:   # too short, not even have {n_stride} exons
                check.append('*')
        else:
            check.append(o)
    return check


def translatability_check(uid,orf):
    check = []
    ensgid = uid.split(':')[0]
    df_certain = df_exonlist.loc[df_exonlist['EnsGID']==ensgid,['EnsPID','Exons']]
    for o,p in zip(orf,df_certain['EnsPID']):
        if o != 'unrecoverable' and o != '':
            result = check_translation(ensgid,p)
            check.append(result)
        else:
            check.append(o)
    return check

def check_translation(EnsGID,EnsPID):
    if EnsPID == 'None':    # usually from RefSeq dataset
        result = '*'
    elif 'PEP' in EnsPID:   # ENSP854949-PEP
        result = '*'
    else:   # protein coding gene or NMD
        pepAnno = dict_biotype[EnsGID]  #{ENSP:anno,ENSP:anno}
        if pepAnno[EnsPID] == 'protein_coding': result = '#'
        else: result = '*'
    return result


