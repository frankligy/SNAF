#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import matplotlib.pyplot as plt
import subprocess
from io import StringIO
from copy import deepcopy,copy




'''
this script is to query the binding affinity of a peptide (9-10)
'''

'''
part I: using netMHCpan4.1b

1. download from https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1
2. tar -xzvf netMHCpan-4.1b.Linux.tar.gz
3. download the data.tar.gz to the netMHCpan-4.1 folder, untar it using the above command, delete the data.tar.gz
4. modify the netMHCpan file, instruction as mentioned in the readme:

      a. At the top of the file  locate the part labelled  "GENERAL SETTINGS:
         CUSTOMIZE TO YOUR SITE"  and set  the 'NMHOME' variable  to the full
	     path to the 'netMHCpan-4.1' directory on your system;

      b. Set TMPDIR to the full path to the temporary directory of you choice. It must
         be user-writable. You may for example set it to $NMHOME/tmp (and create
         the tmp folder in the netMHCpan-4.1 directory).

5. done
'''


def run_netMHCpan(software_path,peptides,hlas,length,cmd_num=1,tmp_dir=None,tmp_name=None):
    # set the default
    if tmp_dir is None:
        tmp_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'scratch')
    if tmp_name is None:
        tmp_name = 'input_{}.pep'.format(os.getpid())
    # reformat/create to the strings that we need
    peptides_path = os.path.join(tmp_dir,tmp_name)
    with open(peptides_path,'w') as f:
        for pep in peptides:
            f.write('{}\n'.format(pep))
    # netMHCpan has a hard limit for 1024 characters for -a, so to be safe, 90 HLA as max for one goal, each input format 11 character 'HLA-B57:01,'
    if len(hlas) <= 90:
        hla_strings = ','.join(hlas)
        if cmd_num == 1:
            reconstruct = ' || '.join(['$3 == ' + '"{}"'.format(pep) for pep in peptides])
            cmd = '{} -p {} -a {} -l {} | awk \'BEGIN {{OFS = "\\t"}} {{if ({}) {{print $3, length($3),$2,$13, $15}}}}\''.format(software_path,peptides_path,hla_strings,length,reconstruct)
            '''
            ../external/netMHCpan-4.1/netMHCpan -p ./test.pep -a HLA-A01:01,HLA-A02:01 -l 9 | awk 'BEGIN {OFS = "\t"} {if ($3 == "AAAWYLWEV" || $3 == "AAGLQDCTM" || $3 == "AARNIVRRA") {print $3, length($3),$2,$13, $15}}'
            '''
        elif cmd_num == 2:
            reconstruct = '\|'.join(peptides)
            cmd = '{} -p {} -a {} -l {} | grep \'{}\' | awk \'BEGIN {{OFS = "\\t"}} {{print $3, length($3),$2,$13, $15}}\''.format(software_path,peptides_path,hla_strings,length,reconstruct)
            '''
            ../external/netMHCpan-4.1/netMHCpan -p ./test.pep -a HLA-A01:01,HLA-A02:01 -l 9 | grep 'AAAWYLWEV\|AAGLQDCTM\|AARNIVRRA' | awk 'BEGIN {OFS = "\t"} {print $3, length($3),$2,$13, $15}'
            '''
        try:
            df = pd.read_csv(StringIO(subprocess.run(cmd,shell=True,capture_output=True).stdout.decode('utf-8')),sep='\t',header=None)
            df.columns = ['peptide','mer','hla','score','identity']
        except:   # no stdout, just no candidates
            df = pd.DataFrame(columns=['peptide','mer','hla','score','identity'])


    else:
        total = len(hlas)   # 91,137,180
        df_store = []
        i = 0 
        while i < total:
            if i + 90 <= total:
                batch_hlas = hlas[i:i+90]
            else:
                batch_hlas = hlas[i:]
            hla_strings = ','.join(batch_hlas)
            if cmd_num == 1:
                reconstruct = ' || '.join(['$3 == ' + '"{}"'.format(pep) for pep in peptides])
                cmd = '{} -p {} -a {} -l {} | awk \'BEGIN {{OFS = "\\t"}} {{if ({}) {{print $3, length($3),$2,$13, $15}}}}\''.format(software_path,peptides_path,hla_strings,length,reconstruct)
                '''
                ../external/netMHCpan-4.1/netMHCpan -p ./test.pep -a HLA-A01:01,HLA-A02:01 -l 9 | awk 'BEGIN {OFS = "\t"} {if ($3 == "AAAWYLWEV" || $3 == "AAGLQDCTM" || $3 == "AARNIVRRA") {print $3, length($3),$2,$13, $15}}'
                '''
            elif cmd_num == 2:
                reconstruct = '\|'.join(peptides)
                cmd = '{} -p {} -a {} -l {} | grep \'{}\' | awk \'BEGIN {{OFS = "\\t"}} {{print $3, length($3),$2,$13, $15}}\''.format(software_path,peptides_path,hla_strings,length,reconstruct)
                '''
                ../external/netMHCpan-4.1/netMHCpan -p ./test.pep -a HLA-A01:01,HLA-A02:01 -l 9 | grep 'AAAWYLWEV\|AAGLQDCTM\|AARNIVRRA' | awk 'BEGIN {OFS = "\t"} {print $3, length($3),$2,$13, $15}'
                '''
            try:
                df = pd.read_csv(StringIO(subprocess.run(cmd,shell=True,capture_output=True).stdout.decode('utf-8')),sep='\t',header=None)
                df.columns = ['peptide','mer','hla','score','identity']
            except:   # no stdout, just no candidates
                df = pd.DataFrame(columns=['peptide','mer','hla','score','identity'])
            df_store.append(df)
            i += 90  
        df = pd.concat(df_store,axis=0)      


    # remove the scratch_pid folder
    os.remove(peptides_path)

    return df



def run_MHCflurry(peptides,hlas):
    peptides = list(peptides)
    tmp_dic_for_alleles= {}
    for index,mhc_ in enumerate(hlas):
        tmp_dic_for_alleles['sample{}'.format(index)] = [mhc_]
    from mhcflurry import Class1PresentationPredictor
    try:
        predictor = Class1PresentationPredictor.load()   # very time consuming
    except:
        cmd = 'mhcflurry-downloads fetch models_class1_presentation'
        subprocess.run(cmd,shell=True)
    result = predictor.predict(peptides=peptides,alleles=tmp_dic_for_alleles,verbose=0) 
    df = result.loc[:,['peptide','best_allele','presentation_percentile']]
    df['mer'] = [len(item) for item in df['peptide']]
    df['identity'] = np.full(shape=df.shape[0],fill_value=None)
    df = df.loc[:,['peptide','mer','best_allele','presentation_percentile','identity']]
    df.rename(columns={'best_allele':'hla','presentation_percentile':'score'},inplace=True)
    return df






