import pandas as pd
import numpy as np
import os,sys
import subprocess
import multiprocessing
import re


def alignment_to_uniprot(df,dict_uni_fa,Ens2ACC,mode):
    col1 = []
    col2 = []
    col3 = []
    for i in range(df.shape[0]):
        # collect all uniprot-curated isoform protein sequence
        target = {}
        EnsID = list(df['UID'])[i].split('|')[0].split(':')[1]
        ACCID = Ens2ACC[EnsID] 
        isoforms = list(dict_uni_fa.keys())  # ['Q9NR97','Q9NR97-2'...]
        for iso in isoforms:        
            if ACCID in iso: 
                seq = dict_uni_fa[iso]
                target[iso] = seq
        # collect all mine-predicted isoform protein sequence
        involve = df['involvement'].tolist()[i]  # [True,False,False,True] indicate which transcript would be a good representative
        match_aa = df['ORFaa'].tolist()[i]
        repre = []
        for idx,aa in enumerate(match_aa):
            if involve[idx] == True:   # only consider aa that is caused by splicing event
                bucket = chop_sequence(aa,10)   # chopping
                subnotes = {}
                for j in range(len(bucket)):   # for each 10mer
                    frag = bucket[j]
                    for key,value in target.items():   # for each curated isoform
                        try: 
                            subnotes[key].append(True) if frag in value else subnotes[key].append(False)
                        except KeyError:
                            subnotes[key] = []
                            subnotes[key].append(True) if frag in value else subnotes[key].append(False)
                for k,m in subnotes.items():   # k is key, m is value
                    if sum(m)==0: subnotes[k].append('notAligned')
                    elif sum(m)==len(m): subnotes[k].append('totallyAligned')
                    else: subnotes[k].append('partiallyAligned')
                repre.append(subnotes)
            elif involve[idx] == False:
                repre.append('Either splicing site is not involved in ORF formation or splicing event doesn\'t occur in this certain transcript')
        col1.append(repre)
        # define what kind of peptide this splicing sites would generate by interogratting each repre list
        identity = []
        for n in repre:
            definition = ''
            if isinstance(n,dict):
                for p in n.values():   #n will be {'P14061':[True,True,False,'partiallyAligned'],'P14061-2':[True,True,False,'partiallyAligned']}
                    if p[-1] == 'totallyAligned':  
                        definition = 'one of already documented isoforms'
                        break
            else: 
                definition = 'Either splicing site is not involved in ORF formation or splicing event doesn\'t occur in this certain transcript'
            if definition == '': definition = 'novel isoform'
            identity.append(definition)
        col2.append(identity)
        # let's see if it is possible to generate any novel isoform


        if mode=='strigent':   # as long as one of possible concatenation results in totallyAligned, then we stop pursuing

            for idx,w in enumerate(identity):
                crystal = True   
                if w == 'one of already documented isoforms': 
                    crystal = False
                    break
            if crystal == True:
                crystal_ = []
                for idx,w in enumerate(identity):
                    if w == 'novel isoform': 
                        query_aa = match_aa[idx]
                        try: result = TMHMM(query_aa,'{0}_{1}_{2}'.format(i,idx,EnsID))  # there's some rare case, it can match up with exising one,
                        # but it doesn't have reasonable ORFs, so the returned ORF prediction is ''
                        except: result = False
                        if result:
                            crys = (True,idx)  # we need look into that
                            crystal_.append(crys)
                if crystal_ == []: col3.append(False) # no need to look into this event anymore
                else: col3.append(crystal_)
            else:
                col3.append(False)
            
        elif mode=='loose':    # consider every possible novel isoform possibilities
            crystal = []
            for idx,w in enumerate(identity):
                if w == 'novel isoform': 
                    query_aa = match_aa[idx]
                    try: result = TMHMM(query_aa,'{0}_{1}_{2}'.format(i,idx,EnsID))  # there's some rare case, it can match up with exising one,
                    # but it doesn't have reasonable ORFs, so the returned ORF prediction is ''
                    except: result = False
                    if result:
                        crys = (True,idx)  # we need look into that
                        crystal.append(crys)
            if crystal == []: col3.append(False) # no need to look into this event anymore
            else: col3.append(crystal)
            
        
    df['alignment'] = col1
    df['identity'] = col2
    try:df['interest'] = col3
    except: 
        print(col3)
        raise Exception('hi')
    return df


def TMHMM(aa,name):
    # download TMHMM linux version, untar it.
    # change shabang of tmhmm and tmhmmformat.pl to the path of perl 5+ you loaded
    # export the path of tmhmm to $PATH, finsh configuration

    # in my use case, save those for convenience
    # perl: /usr/local/perl/5.20.1/bin/perl
    # tmhmm: /data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin
    if not os.path.exists(os.path.join(outFolder,'TMHMM_temp')): os.makedirs(os.path.join(outFolder,'TMHMM_temp'))
    with open(os.path.join(outFolder,'TMHMM_temp','{}.fasta'.format(name)),'w') as f1:
        f1.write('>peptide_{}\n'.format(name))
        f1.write(aa)
    with open(os.path.join(outFolder,'TMHMM_temp','{}.out'.format(name)),'w') as f2:
        subprocess.run(['tmhmm',os.path.join(outFolder,'TMHMM_temp','{}.fasta'.format(name))],stdout = f2)
    with open(os.path.join(outFolder,'TMHMM_temp','{}.out'.format(name)),'r') as f3:
        next(f3)
        punchline = f3.readline().rstrip('\n').split(' ')
        TMn = int(punchline[-1])
    result = True if TMn > 0 else False
    return result