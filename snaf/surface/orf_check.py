import pandas as pd
import numpy as np
import os,sys
import subprocess
import multiprocessing
import re


def get_exon_seq(exon,ensgid,dict_exonCoords):
    attrs = dict_exonCoords[ensgid][exon]
    strand = attrs[1]
    suffer = attrs[4]
    if strand == '+' and not suffer:
        frag = query_from_dict_fa(attrs[2],attrs[3],ensgid,attrs[1])
    elif strand == '+' and suffer:
        frag = query_from_dict_fa(attrs[2],int(attrs[3]-1),ensgid,attrs[1])
    elif strand == '-' and not suffer:
        frag = query_from_dict_fa(attrs[2],attrs[3],ensgid,attrs[1])
    elif strand == '-' and suffer:
        frag = query_from_dict_fa(int(attrs[2]+1),attrs[3],ensgid,attrs[1])
    return frag

def back_traverse_exonlist(exonlist,ensgid,dict_exonCoords,n_stride):
    exonlist = exonlist.split('|')
    exon_iterator = iter(exonlist[::-1])
    current_subexon = next(exon_iterator,'end')  # arrive the last one
    current_exon,current_index = current_subexon[1:].split('.')   # str(15)  # str(1)
    current_exon,current_index = int(current_exon),int(current_index) # int(15)  # int(1)
    current_subexon_seq = get_exon_seq(current_subexon,ensgid,dict_exonCoords)
    n_traversed_exon = 0
    traversed_bases = len(current_subexon_seq)
    while True:
        if n_traversed_exon >= n_stride:
            break
        subexon = next(exon_iterator,'end')
        if subexon == 'end':
            break
        else:
            exon,index = subexon[1:].split('.')
            exon,index = int(exon),int(index)
            subexon_seq = get_exon_seq(subexon,ensgid,dict_exonCoords)
            traversed_bases += len(subexon_seq)
            if exon < current_exon:   # from E15.1 to E14.2
                n_traversed_exon += 1
            else:
                if index < current_index - 1:     # from E14.5 to E14.3
                    n_traversed_exon += 1
    return traversed_bases,n_traversed_exon

def nmd_check(df_exonlist,uid,full_length,orf,dict_exonCoords,n_stride):
    check = []
    ensgid = uid.split(':')[0]
    df_certain = df_exonlist.loc[df_exonlist['EnsGID']=ensgid,['EnsPID','Exons']]
    for o,t,el in zip(orf,full_length,df_certain['Exons']):
        if o != 'unrecoverable' and o != '':
            orf_last_base_pos = t.index(o) + (len(o) - 1)    # think about the start pos is 2, len is 3, means forward 2 (3-1) position
            relative_from_end = len(t) - orf_last_base_pos   # draw a simple graph and think about that.
            traversed_bases, n_traversed_exon = back_traverse_exonlist(el,ensgid,dict_exonCoords,n_stride)
            if n_traversed_exon == n_stride:
                if traversed_bases < relative_from_end:
                    check.append('*')
                else:
                    check.append('#')
            else:   # too short, not even have {n_stride} exons
                check.append('*')
        else:
            check.append(o)




def translatability_check(df_exonlist,uid,orf):
    check = []
    ensgid = uid.split(':')[0]
    df_certain = df_exonlist.loc[df_exonlist['EnsGID']=ensgid,['EnsPID','Exons']]
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


