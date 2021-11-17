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
        Then to run maxQuant:
        cd /path/to/folder/have/raw_and_mqpar
        module load mono/5.20.1
        mono --version
        export PATH=$PATH:​/usr/local/mono/5.20.1/bin
        mono /data/salomonis2/LabFiles/Frank-Li/MaxQuant/bin/MaxQuantCmd.exe /path/to/mqpar.xml
        '''
#######################  Here marks the end of maxQuant configuration