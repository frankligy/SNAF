#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import subprocess
import os,sys


sfa = pd.read_csv('../altanalyze_output/inference/sfa_network.txt',sep='\t',index_col=0)










# ee = pd.read_csv('/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t')
# samples = pd.read_csv('/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/groups.txt',sep='\t',header=None)[0].tolist()
# ee = pd.concat([ee.iloc[:,:11],ee.loc[:,samples]],axis=1)
# ee.to_csv('/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t',index=None)


# clust = pd.read_csv('/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI-clust.txt',sep='\t',index_col=0)
# clust = clust.loc[:,samples]
# clust.to_csv('/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI-clust.txt',sep='\t')

# psi = pd.read_csv('/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI.txt',sep='\t')
# psi = pd.concat([psi.iloc[:,:11],psi.loc[:,samples]],axis=1)
# psi.to_csv('/data/salomonis2/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/test/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI.txt',sep='\t',index=None)
# sys.exit('stop')








# meta_k = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
# file_k = meta_k.index.tolist()

# meta_c = pd.read_csv('../processed_metadata_control.tsv',sep='\t',index_col=0)
# file_c = meta_c.index.tolist()

# file_all = file_k + file_c

# already = []
# for item in pd.read_csv('../md5sum/files_to_check.txt',sep='\t',index_col=0).index:
#     already.append(os.path.basename(item).split('.')[0])

# miss = set(file_all).difference(set(already))

# meta = pd.concat((meta_k,meta_c),axis=0)

# with open('../fastq_miss.txt','w') as f:
#     for item in meta.loc[meta.index.isin(miss),:]['File download URL']:
#         f.write('"{}"\n'.format(item))



























'''
explanation of this block of scratch (2021/10/16):

So the first time download is not completely successful, it should be 2526 fastq to download, but only 2512 get downloaded, so 14 are missing.
Then md5sum check identify 12 are not correct. So I want to first remvoe these 12 wrong fastq on disk, then reconstruct the downlod URL for 12 + 14
= 26 fastq file. Then we can start to download it again.
'''

# comparison = pd.read_csv('../md5sum/comparison.txt',sep='\t')
# remove = comparison.loc[~(comparison['result']),:]['file'].values
# for item in remove:
#     subprocess.run(['rm',item])


# total = []
# with open('../fastq.txt','r') as f:
#     for line in f:
#         total.append(os.path.basename(line.rstrip('"').lstrip('"')).split('.')[0])

# os.chdir('../fastq')
# already = [item.split('.')[0] for item in subprocess.run("for file in *.fastq.gz; do echo $file; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]]

# miss = list(set(total).difference(set(already)))
# meta_k = pd.read_csv('../processed_metadata_KD.tsv',sep='\t',index_col=0)
# meta_c = pd.read_csv('../processed_metadata_control.tsv',sep='\t',index_col=0)
# meta = pd.concat([meta_k,meta_c],axis=0)
# tmp = meta.loc[miss,:]['File download URL'].tolist()
# with open('../fastq_missing.txt','w') as f:
#     for item in tmp:
#         f.write('"{}"\n'.format(item))


'''
End of this block (2021/10/16)
'''
