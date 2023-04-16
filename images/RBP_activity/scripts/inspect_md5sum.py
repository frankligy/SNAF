#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os
import sys


# generate files on which we want to check md5sum
# bash: for file in $(pwd)/*.fastq.gz; do echo $file; done > ../md5sum/files_to_check.txt
# bash: cat files_to_check.txt | xargs -L 1 -P 50 md5sum > md5sum_out.txt

# compare with correct one
output = pd.read_csv('../md5sum/md5sum_out.txt',delim_whitespace=True,header=None)
output.columns = ['md5sum','file']
output['short'] = [os.path.basename(item).split('.')[0] for item in output['file']]

meta1 = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
meta2 = pd.read_csv('../processed_metadata_control.tsv',sep='\t',index_col=0)
meta = pd.concat([meta1,meta2],axis=0)
output['correct'] = output['short'].map(meta['md5sum'].to_dict()).values
output['result'] = [True if output.iloc[i].loc['md5sum'] == output.iloc[i].loc['correct'] else False for i in range(output.shape[0])]
output.to_csv('../md5sum/comparison.txt',sep='\t',index=None)






