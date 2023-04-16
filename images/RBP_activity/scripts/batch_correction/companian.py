#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os

import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from functools import reduce
from copy import deepcopy



# # normalize the junction count, CPM, this step can be done in R, not too slow
# count_file = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/ExpressionInput/counts.original.txt'
# count_norm_file = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/ExpressionInput/counts.original.norm.txt'
# df = pd.read_csv(count_file,sep='\t',index_col=0)
# data = df.values
# sf = data.sum(axis=0).reshape(1,-1) / 1e6
# post = data / sf
# post = pd.DataFrame(data=post,index=df.index,columns=df.columns)
# post.to_csv(count_norm_file,sep='\t')


# # process post-correction
count_norm_post_correct_file = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/ExpressionInput/counts.original.limma.txt'
count_norm_post_correct_processed_file = '/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/ExpressionInput/counts.original.limma.processed.txt'
data = pd.read_csv(count_norm_post_correct_file,sep='\t',index_col=0)
data.index.name = 'AltAnalyze_ID'
data = data.clip(lower=0)
data = data * 25
data.to_csv(count_norm_post_correct_processed_file,sep='\t')






