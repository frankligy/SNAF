#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt

# for biopython, pip install biopython
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq

'''
the pipeline allows seamless integration with all mainstream splicing detection algorithms
'''

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


def construct_dict_exonlist(transcript_db):
    # this will take the mRNA-exonid file, construct to a lookup table
    # ENSG: [E1.1|E3.4|E4.5,E1.1|E3.4|E4.6]
    df_exonlist = pd.read_csv(transcript_db,sep='\t',header=None,names=['EnsGID','EnsTID','EnsPID','Exons'])  # index is number
    dic = df_exonlist.groupby(by='EnsGID')['Exons'].apply(lambda x:x.tolist()).to_dict()
    return dic








