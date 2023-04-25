#!/data/salomonis2/LabFiles/Frank-Li/mhc/test_AltNeo/AltNeo_env/bin/python3.6

# add quote so curl can work
# with open('/Volumes/salomonis2/NCI-R01/Signatures/encode_eclip_bed/raw/K562/bed_url.txt','r') as f1,\
#     open('/Volumes/salomonis2/NCI-R01/Signatures/encode_eclip_bed/raw/K562/bed_url_double_quoted.txt','w') as f2:
#     for line in f1:
#         now = '"' + line.rstrip('\n') + '"' + '\n'
#         f2.write(now)

import pandas as pd
import numpy as np
import subprocess

meta = pd.read_csv('/data/salomonis2/NCI-R01/Signatures/encode_eclip_bed/raw/K562/metadata.tsv',sep='\t')
need = meta.loc[(meta['File type']=='bed')&(meta['File assembly']=='GRCh38')&(meta['Biological replicate(s)']!='1, 2'),:]
col1 = [item + '.bed' for item in need['File accession']]
col2 = [item[:-6] for item in need['Experiment target']]
col3 = [item.split('_')[0] for item in need['Technical replicate(s)']]
col4 = [item[0] + '_' + item[1] for item in zip(col2,col3)]
col5 = [item[0] + '_' + item[1] + '.bed' for item in zip(need['File accession'],col4)]
final = pd.DataFrame({'old':col1,'new':col5,'target':col2})
lookup = final.groupby(by='target')['old'].apply(lambda x:x.tolist())

threshold = 0.1
for i in range(len(lookup)):
    target = lookup.index[i]
    files = lookup.values[i]
    file1 = files[0]
    file2 = files[1]
    subprocess.run(['idr','--samples',file1,file2,'--input-file-type','narrowPeak','--output-file','../idr_bed/{}.bed'.format(target),
                    '--idr-threshold', str(threshold),'--allow-negative-scores'])