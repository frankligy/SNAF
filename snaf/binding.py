#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import matplotlib.pyplot as plt
import subprocess
from io import StringIO

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
    # create the tmp_dir folder if not exist
    if not os.path.exists(tmp_dir):
        os.mkdir(tmp_dir)
    # reformat/create to the strings that we need
    peptides_path = os.path.join(tmp_dir,tmp_name)
    with open(peptides_path,'w') as f:
        for pep in peptides:
            f.write('{}\n'.format(pep))
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
    df = pd.read_csv(StringIO(subprocess.run(cmd,shell=True,capture_output=True).stdout.decode('utf-8')),sep='\t',header=None)
    df.columns = ['peptide','mer','hla','score','identity']
    return df










