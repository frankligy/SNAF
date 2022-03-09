#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt
import xmltodict
import copy
from tqdm import tqdm

from Bio.SeqIO.FastaIO import SimpleFastaParser

########################## Following is for manipulating fasta db files
def chop_normal_pep_db(fasta_path,output_path,mers,allow_duplicates):
    # for human genome in uniprot, 9-10mer, remove duplicates will decrease from 44,741,578 to 41,638,172
    if allow_duplicates:
        with open(fasta_path,'r') as in_handle, open(output_path,'w') as out_handle:
            for title,seq in tqdm(SimpleFastaParser(in_handle)):
                count = 0
                length = len(seq)
                for i in range(length):
                    for mer in mers:
                        if i+mer <= length:
                            out_handle.write('>{}_{}_{}\n'.format(title,mer,count))
                            out_handle.write('{}\n'.format(seq[i:i+mer]))
                            count += 1
    else:
        with open(fasta_path,'r') as in_handle, open(output_path,'w') as out_handle:
            existing = set()
            for title,seq in tqdm(SimpleFastaParser(in_handle)):
                count = 0
                length = len(seq)
                for i in range(length):
                    for mer in mers:
                        if i+mer <= length:
                            subseq = seq[i:i+mer]
                            if subseq not in existing:                                
                                out_handle.write('>{}_{}_{}\n'.format(title,mer,count))
                                out_handle.write('{}\n'.format(subseq))
                                existing.add(subseq)
                                count += 1        

def compare_two_fasta(fa1_path,fa2_path,outdir='.',write_comm=False,write_unique1=False,write_unique2=False,prefix=''):
    seq1,seq2 = set(),set()
    seq1_l,seq2_l = [],[]
    seq1_t,seq2_t = [],[]
    for i,path in enumerate([fa1_path,fa2_path]):
        with open(path,'r') as f:
            for t,s in SimpleFastaParser(f):
                if i == 0:
                    if s not in seq1:
                        seq1.add(s)
                        seq1_l.append(s)
                        seq1_t.append(t)
                elif i == 1:
                    if s not in seq2:
                        seq2.add(s)
                        seq2_l.append(s)
                        seq2_t.append(t)
    comm = seq1.intersection(seq2)
    unique1 = seq1.difference(seq2)
    unique2 = seq2.difference(seq1)
    print('common:{}\nunique1:{}\nunique2:{}'.format(len(comm),len(unique1),len(unique2)))
    if write_comm:
        print('writing comm')
        with open(os.path.join(outdir,'{}comm.fasta'.format(prefix)),'w') as f:
            for item in comm:
                seq1_list = seq1_l
                item_t = seq1_t[seq1_list.index(item)]
                f.write('>{}\n{}\n'.format(item_t,item))
    if write_unique1:
        print('writing unique1')
        with open(os.path.join(outdir,'{}unique1.fasta'.format(prefix)),'w') as f:
            for item in unique1:
                seq1_list = seq1_l
                item_t = seq1_t[seq1_list.index(item)]
                f.write('>{}\n{}\n'.format(item_t,item))
    if write_unique2:
        print('writing unique2')
        with open(os.path.join(outdir,'{}unique2.fasta'.format(prefix)),'w') as f:
            for item in unique2:
                seq2_list = seq2_l
                item_t = seq2_t[seq2_list.index(item)]
                f.write('>{}\n{}\n'.format(item_t,item))




def remove_redundant(fasta_path,out_path):
    existing = set()
    original_count = 0
    with open(fasta_path,'r') as f1, open(out_path,'w') as f2:
        for t,s in SimpleFastaParser(f1):
            original_count += 1
            if s not in existing:
                f2.write('>{}\n{}\n'.format(t,s))
                existing.add(s)
    print('reduce from {} down to {}'.format(original_count,len(existing)))
                

######################## Marks the end of fasta db manipulation

            
#######################  following functions are all for setting mqpar.xml file for maxQuant
def add_database_file(doc,dbs):
    for db in dbs:
        try:
            template = copy.deepcopy(doc['MaxQuantParams']['fastaFiles']['FastaFileInfo'][-1]) # non-first time
        except:
            template = copy.deepcopy(doc['MaxQuantParams']['fastaFiles']['FastaFileInfo']) # first time
        finally:
            template['fastaFilePath'] = db
            template['identifierParseRule'] = '>(.*)'
            template['descriptionParseRule'] = '>(.*)'
            try:
                doc['MaxQuantParams']['fastaFiles']['FastaFileInfo'].append(template)
            except:
                doc['MaxQuantParams']['fastaFiles']['FastaFileInfo'] = []
                doc['MaxQuantParams']['fastaFiles']['FastaFileInfo'].append(template)
    return doc

def numThreads(doc,thread):
    doc['MaxQuantParams']['numThreads'] = str(thread)
    return doc

def add_input_files(doc,inputs):
    experiment = 1
    for inp in inputs:
        # 'filePaths'
        if not type(doc['MaxQuantParams']['filePaths']['string']) == list:
            doc['MaxQuantParams']['filePaths']['string'] = []
            doc['MaxQuantParams']['filePaths']['string'].append(inp)
        else:
            doc['MaxQuantParams']['filePaths']['string'].append(inp)
        
        # 'experiments'
        if not type(doc['MaxQuantParams']['experiments']['string']) == list:
            doc['MaxQuantParams']['experiments']['string'] = []
            doc['MaxQuantParams']['experiments']['string'].append(experiment)
        else:
            doc['MaxQuantParams']['experiments']['string'].append(experiment)
        experiment += 1
        
        # 'fractions'
        if not type(doc['MaxQuantParams']['fractions']['short']) == list:
            doc['MaxQuantParams']['fractions']['short'] = []
            doc['MaxQuantParams']['fractions']['short'].append('32767')
        else:
            doc['MaxQuantParams']['fractions']['short'].append('32767')
            
        # 'ptms'
        if not type(doc['MaxQuantParams']['ptms']['boolean']) == list:
            doc['MaxQuantParams']['ptms']['boolean'] = []
            doc['MaxQuantParams']['ptms']['boolean'].append('False')
        else:
            doc['MaxQuantParams']['ptms']['boolean'].append('False')
            
        # 'paramGroupIndices'
        if not type(doc['MaxQuantParams']['paramGroupIndices']['int']) == list:
            doc['MaxQuantParams']['paramGroupIndices']['int'] = []
            doc['MaxQuantParams']['paramGroupIndices']['int'].append('0')
        else:
            doc['MaxQuantParams']['paramGroupIndices']['int'].append('0')
            
        # 'referenceChannel'
        if not type(doc['MaxQuantParams']['referenceChannel']['string']) == list:
            doc['MaxQuantParams']['referenceChannel']['string'] = []
            doc['MaxQuantParams']['referenceChannel']['string'].append(None)
        else:
            doc['MaxQuantParams']['referenceChannel']['string'].append(None)
    
    return doc

def change_enzymes(doc,enzymes,mode):

    '''
    mode 0: specific
    mode 3: semi-specific
    mode 4: unspecific
    mode 5: no digestion
    '''

    doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymeMode'] = str(mode)

    if enzymes == None:
        doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymes'] = None
    else:
        for enzyme in enzymes:
            if not type(doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymes']['string']) == list:
                doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymes']['string'] = []
                doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymes']['string'].append(enzyme)
            else:
                doc['MaxQuantParams']['parameterGroups']['parameterGroup']['enzymes']['string'].append(enzyme)
     
    return doc           


def change_fdr(doc,protein_fdr,peptide_fdr,site_fdr):
    doc['MaxQuantParams']['proteinFdr'] = protein_fdr
    doc['MaxQuantParams']['peptideFdr'] = peptide_fdr
    doc['MaxQuantParams']['siteFdr'] = site_fdr
    return doc

def change_contaminants(doc,include):
    doc['MaxQuantParams']['includeContaminants'] = str(include)
    return doc

def change_length(doc,minPepLen,maxPeptideMass,minPeptideLengthForUnspecificSearch,maxPeptideLengthForUnspecificSearch):
    doc['MaxQuantParams']['minPepLen'] = minPepLen
    doc['MaxQuantParams']['maxPeptideMass'] = maxPeptideMass
    doc['MaxQuantParams']['minPeptideLengthForUnspecificSearch'] = minPeptideLengthForUnspecificSearch
    doc['MaxQuantParams']['maxPeptideLengthForUnspecificSearch'] = maxPeptideLengthForUnspecificSearch
    return doc


def set_maxquant_configuration(dbs,n_threads,inputs,enzymes,enzyme_mode,outdir,
                               outname='mqpar.xml',protein_fdr=0.01,peptide_fdr=0.01,site_fdr=0.01,include_contaminants=True,
                               minPepLen=9,maxPeptideMass=4600,minPeptideLengthForUnspecificSearch=9,maxPeptideLengthForUnspecificSearch=10):
    with open (os.path.join(os.path.dirname(__file__),'mqpar.xml'),'r') as fr:
        doc = xmltodict.parse(fr.read()) 
        doc = add_database_file(doc,dbs)
        doc = numThreads(doc,n_threads)
        doc = add_input_files(doc,inputs)        
        doc = change_enzymes(doc,enzymes,enzyme_mode)
        doc = change_fdr(doc,protein_fdr=protein_fdr,peptide_fdr=peptide_fdr,site_fdr=site_fdr)
        doc = change_contaminants(doc,include_contaminants)
        doc = change_length(doc,minPepLen,maxPeptideMass,minPeptideLengthForUnspecificSearch,maxPeptideLengthForUnspecificSearch)
        a = xmltodict.unparse(doc,pretty=True,encoding='utf-8')
        a = a.replace('&gt;','>')
        with open(os.path.join(outdir,outname),'w') as fw:
            fw.write(a)

        '''

        Download the maxQuant (working version, 1.6.14 and 2.0.3.1) zip from official website, on windowns, just lanuch exe, follow the instructions,
        .Net framework should be fine, if .Net core 3.1 is needed, just install it and verify it in cmd. On Linux, Using below instructions.

        cd /path/to/folder/have/raw_and_mqpar
        module load mono/5.20.1
        mono --version
        export PATH=$PATH:​/usr/local/mono/5.20.1/bin
        mono /path/to/folder/bin/maxquantCMD.exe /path/to/mqpar.xml   

        '''
#######################  Here marks the end of maxQuant configuration