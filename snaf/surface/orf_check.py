import pandas as pd
import numpy as np
import os,sys
import subprocess
import multiprocessing
import re


def nmd(df_exonlist,ensgid,orf):
    check = []
    df_certain = df_exonlist.loc[df_exonlist['EnsGID']=ensgid,['EnsPID','Exons']]
    for o,el in zip(orf,df_certain['Exons']):
        if o != 'unrecoverable' and o != '':
            seires = build_sorted_exons(ensgid,el)
            num_exon = len(series) - 1
            orf_end_position = 444  # placeholder
            residing = bisect.bisect_left(series,orf_end_pos)   # which exon it resides on
            if residing <= num_exon-2:
                check.append('*')
            else:
                check.append('#')
        else:
            check.append(o)


def build_sorted_exons(EnsGID,exonlists):  # E1.2|E1.3|E2.3|E3.4
    series = []   # store sorted position for each exon
    start_exon = exonlists.split('|')[0]
    strand = dict_exonCoords[EnsGID][start_exon][1]
    if strand == '+':
        start = dict_exonCoords[EnsGID][start_exon][2]
    else:
        start = dict_exonCoords[EnsGID][start_exon][3]  # negative strand, the right most position will be the start, also the largest number

    # final structure {'EnsID':{E1:[chr,strand,start,end],E2:[chr,strand,start,end]}}   
    exonlist = exonlists.split('|')
    dict_judge = {}
    for j in range(len(exonlist)):
        coords = dict_exonCoords[EnsGID][exonlist[j]]
        strand = coords[1]
        judge = check_exonlist_general(exonlist,j,strand)
        dict_judge[exonlist[j]] = judge
    
        
    dic = {}
    for subexon in exonlist:
        exon_num = subexon.split('.')[0]
        subexon_num = subexon.split('.')[1]
        if exon_num in dic:
            dic[exon_num].append(subexon_num)
        else:
            dic[exon_num] = []
            dic[exon_num].append(subexon_num)  # E14 > [1,2,4,5]
    accum = 0
    for exon,sub in dic.items():
        incre,position = check_consecutive(exon,sub,dict_judge,EnsGID,strand,accum)
        accum += incre
        series.extend(position)
    series.sort()   # ascending order [5,9,15,...]
    series = [0]+series
    return series

def check_consecutive(exon,sub,dict_judge,EnsGID,strand,accum):   # E14  > [1,2,4,5]
    #print(exon,sub,dict_judge,accum)
    position = []
    lis_int = [int(x) for x in sub]
    diff1 = np.diff(lis_int,1)   # array([1,2,1])
    diff1 = [int(x)-1 for x in diff1]    # [0,1,0]
    split = np.nonzero(diff1)[0].tolist()  # if pos=1, it means in original list, after index 1 will have a breaking point
    #print(split)
    if split:   # have breaking point
        split = [y + 1 for y in split]    
        # lis_int contains original list, split contains all the indices that identical to the first one in each subgroup
        result=[lis_int[i:j] for i,j in zip([0]+split,split+[None])] 
        for chunk in result:  # chunk[1,2], chunk[4,5]
            query_s = str(exon)+'.'+str(chunk[0])
            query_e = str(exon)+'.'+str(chunk[-1])
            if strand=='+':
                start = dict_exonCoords[EnsGID][query_s][2]
                end = dict_exonCoords[EnsGID][query_e][3] if dict_judge[query_e] else int(dict_exonCoords[EnsGID][query_e][3])-1
                relaPos = int(end) - int(start) + 1   # think 5-1=4, but 5 will be the 5th one
                position.append(relaPos+accum)
            elif strand == '-':
                start = dict_exonCoords[EnsGID][query_s][3]
                end = dict_exonCoords[EnsGID][query_e][2] if dict_judge[query_e] else int(dict_exonCoords[EnsGID][query_e][2])+1
                relaPos = int(start) - int(end) + 1
                position.append(relaPos+accum)
    else:   # E15 > [1,2,3]   3 consecutive 
        query_s = str(exon) + '.' + str(sub[0])
        query_e = str(exon) +'.'+ str(sub[-1])
        #print(query_s,query_e)
        
        if strand=='+':
            start = dict_exonCoords[EnsGID][query_s][2]
            end = dict_exonCoords[EnsGID][query_e][3] if dict_judge[query_e] else int(dict_exonCoords[EnsGID][query_e][3])-1
            relaPos = int(end) - int(start) + 1   # think 5-1=4, but 5 will be the 5th one
            position.append(relaPos+accum)
        elif strand=='-': 
            start = dict_exonCoords[EnsGID][query_s][3]
            end = dict_exonCoords[EnsGID][query_e][2] if dict_judge[query_e] else int(dict_exonCoords[EnsGID][query_e][2])+1
            relaPos = int(start) - int(end) + 1
            position.append(relaPos+accum)
        #print(relaPos)

                
    return relaPos,position


def translatability(df_exonlist,ensgid,orf):
    check = []
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


