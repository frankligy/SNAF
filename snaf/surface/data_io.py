import pandas as pd
import numpy as np
import os,sys
import subprocess
import multiprocessing
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq



def read_uniprot_seq(path):
    dict_fa = {}
    with open(path,'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            ensgid = title.split('|')[3]
            accid = title.split('|')[1]
            try:
                dict_fa[ensgid][accid] = seq
            except KeyError:
                dict_fa[ensgid] = {}
                dict_fa[ensgid][accid] = seq
    return dict_fa   

def fasta_to_dict(path):
    '''
    Let's talk about the fasta file
    >ENSG|chro|start|end
    seq
    the start and end always correspond to xth base in forward strand,
    however, if the ENSG falls into backward strand, the seq it stored is actually the 
    backward strand from its own 5' - 3'.
    '''
    dict_fa = {}  # {ENSID, [chro,start,end,seq]}
    with open(path,'r') as in_handle:
        for title,seq in SimpleFastaParser(in_handle):
            temp_list = []
            EnsID = title.split('|')[0]
            chro = title.split('|')[1]
            start = title.split('|')[2]
            end = title.split('|')[3]
            temp_list=[chro,start,end,seq]
            dict_fa[EnsID] = temp_list
    return dict_fa

def exonCoords_to_dict(path):
    '''
    1. the start and end always forward strand
    2. to clarify the overhang issue, the issue is every middle subexon, its end coord need to be backtracted by 1.
    However, this is different operationally in + and - strand. positive strand is to substract the end by 1 (for middle subexon).
    the negative strand is to add the start by 1 (for middle subexon.)
    '''
    coords=[]
    dict_exonCoords={} # {'EnsID':{E1.1:[chr,strand,start,end,suffer]}} 
    with open(path,'r') as file:
        next(file)
        for line in file:
            items = line.split('\t')
            coords=(items[2],items[3],items[4],items[5],items[10].rstrip('\n'))
            if items[0] in dict_exonCoords:
                dict_exonCoords[items[0]][items[1]] = coords
            else:
                dict_exonCoords[items[0]] = {}
                dict_exonCoords[items[0]][items[1]] = coords
    return dict_exonCoords

def biotype(df):
    # make biotype df to a dictionary
    dic = {} # {EnsGID:{EnsPID:Anno,EnsPID:Anno}}
    for i in range(df.shape[0]):
        EnsGID = df.iat[i,0]
        EnsPID = df.iat[i,1]
        Anno = df.iat[i,2]
        try:
            dic[EnsGID][EnsPID] = Anno
        except KeyError:
            dic[EnsGID] = {EnsPID:Anno}
    return dic
    






