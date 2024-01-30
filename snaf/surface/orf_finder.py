import pandas as pd
import numpy as np
import os,sys
import subprocess
import multiprocessing
import re
from decimal import Decimal as D
from collections import Counter
# for biopython, pip install biopython
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq



def score_GC(sequence):
    length_seq = len(sequence)
    counter = Counter(sequence)
    GC_percent = (counter.get('G',0) + counter.get('C',0)) / length_seq
    return GC_percent

def score_coding_potential(sequence):
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
    min_freq = min(list(usage_dict.values()))
    max_freq = max(list(usage_dict.values()))
    norm_usage_dict = {k:(v-min_freq)/(max_freq-min_freq) for k,v in usage_dict.items()}      
    length_seq = len(sequence)
    num_triplet = length_seq / 3
    i = 0   
    score = 0
    while i + 2 < length_seq:
        triplet = sequence[i:i+3:1]
        score_tri = norm_usage_dict[triplet]
        score += score_tri
        i += 3
    score_potential = score/num_triplet # scale by the number of triplet in the sequence
    return score_potential

def orf2pep(orf):
    assert len(orf) % 3 == 0
    pep = str(Seq(orf).translate(to_stop=False))
    assert '*' not in pep
    return pep

def transcript2orf(cdna):
    candidate_orfs = []
    p_start = re.compile(r'ATG')
    p_end = re.compile(r'(TAA|TGA|TAG)')
    ms = list(re.finditer(p_start,cdna))   
    if len(ms) > 0:
        for s in ms:
            s_pos = int(s.start())
            me = list(re.finditer(p_end,cdna))
            if len(me) > 0:
                for e in me:
                    e_pos = int(e.start())
                    if s_pos < e_pos and s_pos % 3 == e_pos % 3:
                        orf = cdna[s_pos:e_pos]
                        valid = True
                        mo = list(re.finditer(p_end,orf))   # still need to check whether there is stop codon in the orf
                        if len(mo) > 0:
                            for o in mo:
                                o_pos = int(o.start())
                                if o_pos % 3 == 0:
                                    valid = not valid
                                    break
                        if valid:
                            candidate_orfs.append(orf)
            else:  # maybe the stop codon not in the cdna, the last reading codon is the last three bases
                orf = cdna[s_pos:]
                valid = True
                if len(orf) % 3 != 0:
                    valid = not valid
                mo = list(re.finditer(p_end,orf))
                if len(mo) > 0:
                    for o in mo:
                        o_pos = int(o.start())
                        if o_pos % 3 == 0:
                            valid = not valid
                            break
                if valid:
                    candidate_orfs.append(orf)
    return candidate_orfs



def prioritize_orf(candidate_orfs,min_len=30*3,tol_len=8*3):
    max_orf = ''
    max_length = 0
    max_score = 0
    for orf in candidate_orfs:
        if len(orf) < min_len:
            continue
        score = score_GC(orf) + score_coding_potential(orf)
        if len(orf) - max_length > tol_len:
            max_length = len(orf)
            max_score = score
            max_orf = orf
        else:
            if len(orf) > max_length and score > max_score:
                max_length = len(orf)
                max_score = score
                max_orf = orf  
    return max_orf            


