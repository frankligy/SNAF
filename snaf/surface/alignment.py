import pandas as pd
import numpy as np
import os,sys
import subprocess
import multiprocessing
import re

def alignment_to_uniprot(orf,uid,dict_uni_fa,tmhmm=False,software_path=None):
    ensgid = uid.split(':')[0]
    isoforms = dict_uni_fa[ensgid]  # {acc1:seq,acc1-2:seq}
    results = []
    for o in orf:
        if o != 'unrecoverable' and o != '':
            final_notes = []
            for iso,iso_seq in isoforms.items():
                if len(iso_seq) != len(o):
                    conclusion = True
                else:
                    notes = []                    
                    for i, mer in enumerate(chop_sequence(o,10)):
                        if mer in iso_seq:
                            notes.append(True)  # here True means align
                        else:
                            notes.append(False)
                    conclusion = not all(notes)   # here True means not existing before
                final_notes.append(conclusion)
            decision = all(final_notes)  # here True means this isoform is novel
            if tmhmm:
                is_cross_membrane = TMHMM(o,software_path=software_path)
                decision = all([decision,is_cross_membrane])
            results.append(decision)
        else:
            results.append(o)
    return results
    

def chop_sequence(seq,kmer):   # how to splice sequence, elegant way to use range
    frag_bucket = []
    for i in range(0,len(seq),kmer):
        if i + kmer <= len(seq):
            frag_bucket.append(seq[i:i+kmer])
        else:
            frag_bucket.append(seq[i:])
    return frag_bucket


def TMHMM(aa,software_path):
    # download TMHMM linux version, untar it.
    # change shabang of tmhmm and tmhmmformat.pl to the path of perl 5+ you loaded
    # export the path of tmhmm to $PATH, finsh configuration

    # in my use case, save those for convenience
    # perl: /usr/local/perl/5.20.1/bin/perl   # no need to load perl when running tmhmm
    # tmhmm: /data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin
    name = os.getpid()
    int_file_path = os.path.join(os.path.dirname(os.path.dirname(__file__)),'scratch','{}.fasta'.format(name))
    with open(int_file_path,'w') as f1:
        f1.write('>peptide_{}\n'.format(name))
        f1.write(aa)
    '''
    # sp|Q9NR97|TLR8_HUMAN Length: 942
    # sp|Q9NR97|TLR8_HUMAN Number of predicted TMHs:  1
    # sp|Q9NR97|TLR8_HUMAN Exp number of AAs in TMHs: 23.20959
    # sp|Q9NR97|TLR8_HUMAN Exp number, first 60 AAs:  0.03874
    # sp|Q9NR97|TLR8_HUMAN Total prob of N-in:        0.00259
    sp|Q9NR97|TLR8_HUMAN	TMHMM2.0	outside	     1   825
    sp|Q9NR97|TLR8_HUMAN	TMHMM2.0	TMhelix	   826   848
    sp|Q9NR97|TLR8_HUMAN	TMHMM2.0	inside	   849   942
    '''
    lines = subprocess.run('{} {}'.format(software_path,int_file_path),shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    n_cross = int(lines[1].split('  ')[-1])
    result = True if n_cross > 0 else False
    return result