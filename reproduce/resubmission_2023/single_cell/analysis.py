#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
sys.path.insert(0,'/data/salomonis2/software')
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import anndata as ad
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import subprocess

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# # preprocess
# df = pd.read_csv('GSE72056_melanoma_single_cell_revised_v2.txt',sep='\t',header=[0,1,2,3],index_col=0)
# '''
# 4 levels of columns 
# Cell
# tumor
# malignant(1=no,2=yes,0=unresolved)
# non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)
# '''
# all_cells = df.columns.get_level_values('Cell').tolist()
# col = []
# for item in all_cells:
#     if 'CY94_CD45NEG_CD90POS_2' in item:  # CY94_CD45NEG_CD90POS_2_C02_S26_comb to CY94-CD45NEG-CD90POS-2-ALL-C02_S26_comb
#         former = '-'.join(item.split('_')[:-2])
#         latter = '_'.join(item.split('_')[-2:])
#         new = '_'.join([former,latter]) # CY94-CD45NEG-CD90POS-2-C02_S26_comb 
#         new = '-'.join(new.split('-')[:-1]) + '-ALL-' + new.split('-')[-1]
#     else:
#         lis = item.split('_')
#         if len(lis) > 2:  # cy79-p1-CD45-pos-PD1-neg-AS-C1-R2-D07-S523-comb to cy79-p1-CD45-pos-PD1-neg-AS-C1-R2-D07_S523_comb
#             former = '-'.join(item.split('_')[:-2])
#             latter = '_'.join(item.split('_')[-2:])
#             new = '_'.join([former,latter])
#         elif len(lis) == 2:  # Cy67-CD45pos-S2-A10_S10 to Cy67-CD45pos-S2-A10_S10_L001
#             new = item + '_L001'
#         elif len(lis) == 1:   # Cy72_CD45_H02_S758_comb to Cy72-CD45-H02_S758_comb
#             former = '-'.join(item.split('-')[:-2])
#             latter = '_'.join(item.split('-')[-2:])
#             new = '_'.join([former,latter])
#         else:
#             new = item
#     col.append(new)
# all_fastqs = pd.read_csv('all_fastqs.txt',sep='\t',header=None,index_col=0).index.tolist()
# all_fastqs = [item.split('_R1_001')[0] for item in all_fastqs]
# common = list(set(col).intersection(set(all_fastqs)))
# not_map = list(set(col).difference(set(all_fastqs)))
# cond = [True if item in common else False for item in col]
# # with open('not_be_able_to_map.txt','w') as f:
# #     for item in not_map:
# #         f.write('{}\n'.format(item))
# lookup = pd.DataFrame(data={'barcode':all_cells,'fastq':col,'cond':cond})
# lookup.to_csv('lookup.txt',sep='\t')

# # move those mappable fastq to a folder for following analysis
# lookup = pd.read_csv('lookup.txt',sep='\t',index_col=0)
# lookup = lookup.loc[lookup['cond'],:]  # 4456 + 189 = 4645
# for item in tqdm(lookup['fastq'],total=lookup.shape[0]):
#     r1 = item + '_R1_001.fastq.gz'
#     r2 = item + '_R2_001.fastq.gz'
#     cmd1 = 'cp ../aviv_melanoma/{} ./fastqs'.format(r1)
#     cmd2 = 'cp ../aviv_melanoma/{} ./fastqs'.format(r2)
#     subprocess.run(cmd1,shell=True)
#     subprocess.run(cmd2,shell=True)

# # generate all_samples
# lookup = pd.read_csv('lookup.txt',sep='\t',index_col=0)
# lookup = lookup.loc[lookup['cond'],:]
# with open('all_samples.txt','w') as f:
#     for item in lookup['fastq']:
#         f.write('{}\n'.format(item))

# get cell type assignment
df = pd.read_csv('GSE72056_melanoma_single_cell_revised_v2.txt',sep='\t',header=[0,1,2,3],index_col=0)
identity = df.columns.to_frame(index=False).set_index(keys='Cell')
# tumor 
# malignant(1=no,2=yes,0=unresolved) 
# non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)
ct = []
dic = {'1':'Tcell','2':'Bcell','3':'Macrophage','4':'Endothelial','5':'CAF','6':'NK','0':'unknown'}
for row in identity.iterrows():
    row = row[1]
    if row['malignant(1=no,2=yes,0=unresolved)'] == '2':
        ct.append('tumor')
    else:
        ct.append(dic[row['non-malignant cell type (1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)']])
identity['ct'] = ct
identity.to_csv('celltype.txt',sep='\t')

# df.columns = df.columns.get_level_values('Cell').tolist()
# df.to_csv('normed_expression_matrix.txt',sep='\t')


'''hla typing, just choose 19 bam files for that'''
'''seems doesn't work
Maximal read length of 0 bps exceeded. Please remove "#define RAZERS_MEMOPT" in razers.cpp and recompile.
Failed to load reads
Exiting ...
'''
# tmp = identity.drop_duplicates(subset='tumor')
lookup = pd.read_csv('lookup.txt',sep='\t',index_col=0).set_index(keys='barcode')
lookup = lookup['fastq'].to_dict()
# col = [lookup[item] for item in tmp.index]
# tmp.index = col
# tmp = tmp.loc[~tmp['tumor'].isin(['59','65']),:]
# tmp.to_csv('represent_bam_for_typing.txt',sep='\t')

'''hla typing, let's randomly sample 50 cells for each tumor, concat the fastq and see'''
# lookup = pd.read_csv('lookup.txt',sep='\t',index_col=0).set_index(keys='barcode')
# lookup = lookup['fastq'].to_dict()
# for t,sub_df in identity.groupby(by='tumor'):
#     if t == '59' or t == '65':
#         continue
#     else:
#         chosen = sub_df.sample(n=50,replace=False).index
#         chosen = [lookup[item] for item in chosen]
#         os.mkdir('typing_fastq/{}'.format(t))
#         for sample in chosen:
#             subprocess.run(['cp','fastqs/{}_R1_001.fastq.gz'.format(sample),'typing_fastq/{}'.format(t)])
#             subprocess.run(['cp','fastqs/{}_R2_001.fastq.gz'.format(sample),'typing_fastq/{}'.format(t)])


'''just analyze count matrix, separated by ct'''
tumor_cells = [item + '_secondAligned.sortedByCoord.out.bed' for item in identity.loc[identity['ct']=='tumor',:].index] 
t_cells = [item + '_secondAligned.sortedByCoord.out.bed' for item in identity.loc[identity['ct']=='Tcell',:].index] 
b_cells = [item + '_secondAligned.sortedByCoord.out.bed' for item in identity.loc[identity['ct']=='Bcell',:].index] 
macrophage_cells = [item + '_secondAligned.sortedByCoord.out.bed' for item in identity.loc[identity['ct']=='Macrophage',:].index] 
caf_cells = [item + '_secondAligned.sortedByCoord.out.bed' for item in identity.loc[identity['ct']=='CAF',:].index] 
nk_cells = [item + '_secondAligned.sortedByCoord.out.bed' for item in identity.loc[identity['ct']=='NK',:].index] 
endo_cells = [item + '_secondAligned.sortedByCoord.out.bed' for item in identity.loc[identity['ct']=='Endothelial',:].index] 

valid_samples = pd.read_csv('altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t').columns.tolist()[11:]

# # splicing events that are high in tumor versus other cell type
# # tumor versus other
# with open('altanalyze_output/DAS/groups.tumor_other.txt','w') as f:
#     for cell in tumor_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,1,'tumor'))
#     for cell in t_cells + b_cells + macrophage_cells + caf_cells + nk_cells + endo_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,2,'other'))
# with open('altanalyze_output/DAS/comps.tumor_other.txt','w') as f:
#     f.write('1\t2\n')

# # tcell versus other
# with open('altanalyze_output/DAS/groups.tcell_other.txt','w') as f:
#     for cell in t_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,1,'tcell'))
#     for cell in tumor_cells + b_cells + macrophage_cells + caf_cells + nk_cells + endo_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,2,'other'))
# with open('altanalyze_output/DAS/comps.tcell_other.txt','w') as f:
#     f.write('1\t2\n')

# # bcell versus other
# with open('altanalyze_output/DAS/groups.bcell_other.txt','w') as f:
#     for cell in b_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,1,'bcell'))
#     for cell in tumor_cells + t_cells + macrophage_cells + caf_cells + nk_cells + endo_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,2,'other'))
# with open('altanalyze_output/DAS/comps.bcell_other.txt','w') as f:
#     f.write('1\t2\n')

# # macrophage versus other
# with open('altanalyze_output/DAS/groups.macrophage_other.txt','w') as f:
#     for cell in macrophage_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,1,'macrophage'))
#     for cell in tumor_cells + t_cells + b_cells + caf_cells + nk_cells + endo_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,2,'other'))
# with open('altanalyze_output/DAS/comps.macrophage_other.txt','w') as f:
#     f.write('1\t2\n')

# # caf versus other
# with open('altanalyze_output/DAS/groups.caf_other.txt','w') as f:
#     for cell in caf_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,1,'caf'))
#     for cell in tumor_cells + t_cells + b_cells + macrophage_cells + nk_cells + endo_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,2,'other'))
# with open('altanalyze_output/DAS/comps.caf_other.txt','w') as f:
#     f.write('1\t2\n')

# # nk versus other
# with open('altanalyze_output/DAS/groups.nk_other.txt','w') as f:
#     for cell in nk_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,1,'nk'))
#     for cell in tumor_cells + t_cells + b_cells + macrophage_cells + caf_cells + endo_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,2,'other'))
# with open('altanalyze_output/DAS/comps.nk_other.txt','w') as f:
#     f.write('1\t2\n')

# # endo versus other
# with open('altanalyze_output/DAS/groups.endo_other.txt','w') as f:
#     for cell in endo_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,1,'endo'))
#     for cell in tumor_cells + t_cells + b_cells + macrophage_cells + caf_cells + nk_cells:
#         cell = cell.split('_secondAligned.sortedByCoord.out.bed')[0]
#         cell = lookup[cell]
#         cell = cell + '_secondAligned.sortedByCoord.out.bed'
#         if cell in valid_samples:
#             f.write('{}\t{}\t{}\n'.format(cell,2,'other'))
# with open('altanalyze_output/DAS/comps.endo_other.txt','w') as f:
#     f.write('1\t2\n')

# compare with TCGA melanoma junction
jcmq = snaf.JunctionCountMatrixQuery.deserialize('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/after_prediction.p')
ts_junctions = jcmq.valid # 16799
ee = pd.read_csv('altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',sep='\t')
uids = list(set([':'.join(item.split('|')[0].split(':')[1:]) for item in ee['UID']]))  # 34245
common = list(set(ts_junctions).intersection(set(uids)))  # 542

n = []
o = []
for ct in ['tumor','tcell','bcell','macrophage','nk','caf','endo']:
    psi = pd.read_csv('altanalyze_output/AltResults/AlternativeOutput/Events-dPSI_0.0_rawp/PSI.{}_vs_other.txt'.format(ct),sep='\t',index_col=0)
    psi = psi.loc[(psi['dPSI']>0.1)&(psi['adjp']<0.05),:]
    n.append(psi.shape[0])
    uids = list(set([':'.join(item.split('|')[0].split(':')[1:]) for item in psi.index]))
    overlap = list(set(common).intersection(set(uids)))
    o.append(len(overlap)) # no overlap
df = pd.DataFrame(data={'ct':['tumor','tcell','bcell','macrophage','nk','caf','endo'],'number':n,'overlap':o})
df.sort_values(by='number',inplace=True,ascending=False)
df.plot.bar(x='ct',y='number')
plt.savefig('das_by_ct.pdf',bbox_inches='tight')
plt.close()




