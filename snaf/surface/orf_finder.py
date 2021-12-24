import pandas as pd
import numpy as np
import os,sys
import subprocess
import multiprocessing
import re
from decimal import Decimal as D



def rescue_position(pos,manner):
    for m in re.finditer(r'ATG',manner):
        if m.start() % 3 ==0:
            span = len(manner) - m.start()
            protruding = span % 3
            end = -1 - protruding + 1
            frag = manner[m.start():end]
            pos.append(frag)
    return pos   # pos is actually a list of peptides

def pos_to_frags(pos,sequence):
    frag_array = []
    if pos:       
        frag_array.append(sequence[0:pos[0]])
        i = 0
        while i < len(pos)-1:
            frag_array.append(sequence[pos[i]+3:pos[i+1]])
            i += 1
        last_seq = sequence[pos[-1]+3:]
    return frag_array, last_seq  


def score_GC(sequence):
    GC_content = 0
    length_seq = len(sequence)
    for nt in sequence:
        if nt == 'G' or nt == 'C':
            GC_content += 1
    try:
        GC_percent = GC_content / length_seq
    except:
        print('here it the seq:',sequence)
        raise Exception
    return GC_percent

def score_coding_bias(sequence):
    # coding frequency table is from GenScript webpage
    usage_dict = {'TTT':16.9,'TTC':20.4,'TTA':7.2,'TTG':12.6,'TAT':12.0,'TAC':15.6,'TAA':0.7,'TAG':0.5,
                  'CTT':12.8,'CTC':19.4,'CTA':6.9,'CTG':40.3,'CAT':10.4,'CAC':14.9,'CAA':11.8,'CAG':34.6,
                  'ATT':15.7,'ATC':21.4,'ATA':7.1,'ATG':22.3,'AAT':16.7,'AAC':19.5,'AAA':24.0,'AAG':32.9,
                  'GTT':10.9,'GTC':14.6,'GTA':7.0,'GTG':28.9,'GAT':22.3,'GAC':26.0,'GAA':29.0,'GAG':40.8,
                  'TCT':14.6,'TCC':17.4,'TCA':11.7,'TCG':4.5,'TGT':9.9,'TGC':12.2,'TGA':1.3,'TGG':12.8,
                  'CCT':17.3,'CCC':20.0,'CCA':16.7,'CCG':7.0,'CGT':4.7,'CGC':10.9,'CGA':6.3,'CGG':11.9,
                  'ACT':12.8,'ACC':19.2,'ACA':14.8,'ACG':6.2,'AGT':11.9,'AGC':19.4,'AGA':11.5,'AGG':11.4,
                  'GCT':18.6,'GCC':28.5,'GCA':16.0,'GCG':7.6,'GGT':10.8,'GGC':22.8,'GGA':16.3,'GGG':16.4} 
    # do a normaliztion for each triplet, then for all the triplet's sum, divided by the number of triplet
    min_freq = 4.5
    max_freq = 40.8
    norm_usage_dict = {}
    for codon,freq in usage_dict.items():
        norm_usage_dict[codon] = float((D(freq) - D(min_freq)) / (D(max_freq) - D(min_freq)))        
    length_seq = len(sequence)
    num_triplet = length_seq/3
    i = 0   
    score = 0
    while i < length_seq - 2:
        triplet = sequence[i:i+3:1]
        score_tri = norm_usage_dict[triplet]
        score += score_tri
        i += 3
    score_bias = score/num_triplet # scale by the number of triplet in the sequence
    return score_bias 


def transcript2peptide(cdna_sequence):   # actually to ORF
    reading_manners = []
    reading_manners.append(cdna_sequence[0:])
    reading_manners.append(cdna_sequence[1:])
    reading_manners.append(cdna_sequence[2:])
    frag_comp_array = []
    for manner in reading_manners:       
        pos = []
        for m in re.finditer(r'(TAA|TGA|TAG)',manner):   # for multiple instances
            if m.start() % 3 == 0:
                pos.append(m.start())
        if pos == []:
            frag_comp_array.extend(rescue_position(pos,manner))
        else:
            frag_array,last_seq = pos_to_frags(pos,manner)
            for frag in frag_array:
                if 'ATG' not in frag or len(frag) == 0:
                    continue
                else:
                    for n in re.finditer(r'ATG',frag):
                        if (len(frag) - n.start()) % 3 == 0:
                            frag_comp = frag[n.start():]
                            frag_comp_array.append(frag_comp)
                            break   # might have multiple 'ATG' so it is necessary to break when find first 'ATG'
                        else:
                            continue
        # process last_seq:
            for n in re.finditer(r'ATG',last_seq):
                if n.start() % 3 == 0:
                    last_frag = last_seq[n.start():]
                    protruding = len(last_frag) % 3
                    end = -1 - protruding + 1   # python end exclusive, so + 1
                    last_frag_real = last_frag[:end]
                    frag_comp_array.append(last_frag_real)
    #######################  # We think if you only has longer length(0-7) but add_score is not higher than original one, you are FAlSE
    max_seq = ''
    max_length = 0
    max_item_score = 0
    for item in frag_comp_array:
        temp1 = len(item)
        if temp1==0: continue
        else:
            add_score = score_GC(item) + score_coding_bias(item)
            if (temp1 - max_length) >= 8:
                max_length = temp1
                max_item_score = add_score
                max_seq = item
            elif (temp1 - max_length) >= 0 and (temp1 - max_length) < 8:
                if add_score >= max_item_score:
                    max_length = temp1
                    max_item_score = add_score
                    max_seq = item
#           else:
#                print('equal length but less likely to be a true ORF or longer length but less likely to be a true ORF',add_score,max_item_score) 
    max_seq_tran = max_seq
    return max_seq_tran


