#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import os,sys,subprocess
import numpy as np
import pandas as pd


# # venn for functional
# base_dir = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/das_within'
# for folder in ['shRNA_HepG2','shRNA_K562','CRISPR_HepG2','CRISPR_K562']:
#     whole_dir = os.path.join(base_dir,folder)
#     os.chdir(whole_dir)
#     all_sfs = subprocess.run('for file in *.txt; do echo $(basename $file .txt); done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     with open('../../scratch/{}_all_sfs.txt'.format(folder),'w') as f:
#         for sf in all_sfs:
#             f.write('{}\n'.format(sf))

# venn for binding
base_dir = '/data/salomonis2/NCI-R01/Signatures/encode_eclip_bed/raw'
for folder in ['K562','HepG2']:
    whole_dir = os.path.join(base_dir,folder,'idr_bed')
    os.chdir(whole_dir)
    all_sfs = subprocess.run('for file in *.bed; do echo $(basename $file .bed); done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    with open('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/scratch/{}_binding_all_sfs.txt'.format(folder),'w') as f:
        for sf in all_sfs:
            f.write('{}\n'.format(sf))
