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


# get reduced junction
df = snaf.get_reduced_junction_matrix(pc='counts.TCGA-SKCM.txt',pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')

# run SNAF
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
surface.initialize(db_dir=db_dir)

for uid in ['ENSG00000105976:E3.1-E4.1','ENSG00000198053:E7.2-E13.1_1915159','ENSG00000164175:E3.2_33963931-E4.2','ENSG00000092421:E22.1-E24.1_116468915','ENSG00000057019:E5.1-I5.1_98867438','ENSG00000152558:I7.1_102419074-E8.1']:
    sigma, df = snaf.tumor_specificity(uid=uid,method='mle',return_df=True)
    fig,ax = plt.subplots()
    sns.histplot(df['value_cpm'],bins=100,kde=True,ax=ax,stat='density')
    from scipy.stats import halfnorm
    r = halfnorm.rvs(loc=0,scale=sigma,size=2000)
    sns.kdeplot(r,ax=ax,color='orange',clip=(-0.01,5))
    plt.savefig('{}.pdf'.format(uid),bbox_inches='tight')
    plt.close()
sys.exit('stop')

# df = pd.read_csv('test.txt',sep='\t',index_col=0)
# new_df = snaf.add_tumor_specificity_frequency_table(df,method='mean',remove_quote=True,cores=None)
# new_df.to_csv('test1.txt',sep='\t')

# for uid in ['ENSG00000241343:E2.6-E2.9_101392052','ENSG00000175482:E2.22_67352741-E3.1']:
#     snaf.gtex_visual_combine_plotly(uid=uid,outdir='.',norm=False,tumor=df)
#     snaf.gtex_visual_combine_plotly(uid=uid,outdir='.',norm=True,tumor=df)     

# for uid in ['ENSG00000084764:E19.1-E20.1','ENSG00000243290:E1.1-E1.1_89040710','ENSG00000165914:E6.3-E8.1']:
#     snaf.gtex_visual_combine_plotly(uid=uid,outdir='.',norm=False,tumor=df)
#     snaf.gtex_visual_combine_plotly(uid=uid,outdir='.',norm=True,tumor=df)   

sys.exit('stop')
# for uid in ['ENSG00000105976:E3.1-E4.1','ENSG00000198053:E7.2-E13.1_1915159','ENSG00000164175:E3.2_33963931-E4.2','ENSG00000092421:E22.1-E24.1_116468915','ENSG00000057019:E5.1-I5.1_98867438','ENSG00000152558:I7.1_102419074-E8.1']:
#     snaf.gtex_visual_combine_plotly(uid=uid,outdir='.',norm=False,tumor=df)
#     snaf.gtex_visual_combine_plotly(uid=uid,outdir='.',norm=True,tumor=df)

'''T cell neoantigen'''
jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=30,add_control=add_control,outdir='result',filter_mode='maxmin')
sample_to_hla = pd.read_csv('sample_hla.txt',sep='\t',index_col=0)['hla'].to_dict()
hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
jcmq.run(hlas=hlas,outdir='./result')
snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')

sys.exit('stop')

'''B cell neoantigen'''
# membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df)
# surface.run(membrane_tuples,prediction_mode='short_read',outdir='result',tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm')
# surface.run(uids=membrane_tuples,outdir='.',prediction_mode='long_read',gtf='2021UHRRIsoSeq_SQANTI3_filtered.gtf',tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm')
# surface.generate_results(pickle_path='./result/surface_antigen_lr.p',outdir='result',strigency=3,gtf=None,long_read=True)   # './SQANTI-all/collapse_isoforms_classification.filtered_lite.gtf'
# surface.run_dash_B_antigen(pkl='result/surface_antigen.p',prediction_mode='short_read',candidates='result/candidates_3.txt',python_executable='/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7')
# surface.run_dash_B_antigen(pkl='result/surface_antigen_lr.p',prediction_mode='long_read',candidates='result/candidates_3_lr.txt',python_executable='/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7')


# print(snaf.uid_to_coord('ENSG00000164175:E3.2_33963931-E4.2'))
# snaf.gtex_visual_combine('ENSG00000057019:E5.1-I5.1_98867438',norm=False,outdir='result',tumor=df)   # ENSG00000198053:E7.2-E13.1_1915159 ENSG00000164175:E3.2_33963931-E4.2 ENSG00000092421:E22.1-E24.1_116468915 ENSG00000152558:I7.1_102419074-E8.1 ENSG00000105976:E3.1-E4.1 ENSG00000057019:E5.1-I5.1_98867438


'''
output B antigen results for paper
'''
# surface.report_candidates('result/surface_antigen.p','result/candidates_3.txt','result/validation_3.txt','result/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt','short_read','result','sr_str3_report.txt')
# surface.report_candidates('result/surface_antigen.p','result/candidates_4.txt','result/validation_4.txt','result/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt','short_read','result','sr_str4_report.txt')
# surface.report_candidates('result/surface_antigen.p','result/candidates_5.txt','result/validation_5.txt','result/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt','short_read','result','sr_str5_report.txt')
# surface.report_candidates('result/surface_antigen_lr.p','result/candidates_3_lr.txt','result/validation_3_lr.txt','result/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt','long_read','result','lr_str3_report.txt')












### MS validation
# df = pd.read_csv('result/frequency_stage3_verbosity1_uid_gene_symbol.txt',sep='\t',index_col=0)
# df_common = df.loc[df['n_sample']>282,:]
# df_unique = df.loc[df['n_sample']==1,:]
# with open('result/MS_validation_common.fasta','w') as f:
#     for row in df_common.itertuples():
#         peptide, uid = row.Index.split(',')
#         f.write('>{}\n{}\n'.format(uid,peptide))
# with open('result/MS_validation_unique.fasta','w') as f:
#     for row in df_unique.itertuples():
#         peptide, uid = row.Index.split(',')
#         f.write('>{}\n{}\n'.format(uid,peptide))
# for f in ['MS_validation_common','MS_validation_unique']:
#     snaf.proteomics.remove_redundant('./result/{}.fasta'.format(f),'./result/{}_unique.fasta'.format(f))
#     snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot_9_10_mers_unique.fasta',
#                                       fa2_path='./result/{}_unique.fasta'.format(f),outdir='./result',write_unique2=True,prefix='{}_'.format(f))

db_common = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/fasta/MS_validation_common_unique2.fasta']
db_unique = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/fasta/MS_validation_unique_unique2.fasta']

'''
common
'''
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# for p in all_patients:
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common/{}'.format(p))
#     all_samples = subprocess.run('for file in *.raw; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     inputs = []
#     pwd = os.getcwd()
#     for s in all_samples:
#         inputs.append(os.path.join(pwd,s))
#     snaf.proteomics.set_maxquant_configuration(dbs=db_common,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                                outdir=pwd,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)


'''
unique
'''
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_unique')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# for p in all_patients:
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_unique/{}'.format(p))
#     all_samples = subprocess.run('for file in *.raw; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     inputs = []
#     pwd = os.getcwd()
#     for s in all_samples:
#         inputs.append(os.path.join(pwd,s))
#     snaf.proteomics.set_maxquant_configuration(dbs=db_unique,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                                outdir=pwd,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)

'''
concat
'''
# with open('result/MS_validation_common_unique2.fasta','r') as f1, open('result/MS_validation_unique_unique2.fasta','r') as f2, open('result/MS_validation_concat_unique2.fasta','w') as f3:
#     # first process common neoantigen fasta
#     for line in f1:
#         if line.startswith('>'):
#             uid = line.rstrip('\n')
#             iden = 'common_neoantigen'
#             new_line = '>{},{}\n'.format(uid,iden)
#             f3.write(new_line)
#         else:
#             f3.write(line)
#     # then process unique neoantigen fasta
#     for line in f2:
#         if line.startswith('>'):
#             uid = line.rstrip('\n')
#             iden = 'unique_neoantigen'
#             new_line = '>{},{}\n'.format(uid,iden)
#             f3.write(new_line)
#         else:
#             f3.write(line)

# db_concat = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/fasta/MS_validation_concat_unique2.fasta']
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# for p in all_patients:
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat/{}'.format(p))
#     all_samples = subprocess.run('for file in *.raw; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     inputs = []
#     pwd = os.getcwd()
#     for s in all_samples:
#         inputs.append(os.path.join(pwd,s))
#     snaf.proteomics.set_maxquant_configuration(dbs=db_concat,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
#                                                outdir=pwd,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)


'''
plot separate
'''
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# fig,ax = plt.subplots()
# n_common = []
# n_unique = []
# from scipy.stats import mannwhitneyu,ks_2samp,ttest_ind,ttest_rel
# all_patients.pop(all_patients.index('Mel-16'))
# for p in all_patients:
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common/{}/combined/txt'.format(p))
#     pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     n = int(pep.shape[0])/114
#     n_common.append(n)
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_unique/{}/combined/txt'.format(p))
#     pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     n = int(pep.shape[0])/25188
#     n_unique.append(n)  
# stats = ttest_ind(n_common,n_unique)
# print(stats)
# number = len(all_patients)
# ax.bar(x=[i for i in range(1,1+3*(number-1)+1,3)],height=n_common,label='common neoantigen')
# ax.bar(x=[i for i in range(2,2+3*(number-1)+1,3)],height=n_unique,label='unique neoantigen')
# xticks = []
# for ic,iu in zip([i for i in range(1,1+3*(number-1)+1,3)],[i for i in range(2,2+3*(number-1)+1,3)]):
#     xticks.append((ic+iu)/2)
# ax.set_xticks(xticks)
# ax.set_xticklabels(all_patients,fontsize=4,rotation=90)
# ax.set_xlabel('patients',fontsize=6)
# ax.set_ylabel('MS recovery rate',fontsize=6)
# ax.legend(loc='upper left',bbox_to_anchor=(1,1),frameon=False)
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# plt.savefig('MS_plot.pdf',bbox_inches='tight')
# plt.close()

'''
plot concat
'''
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# fig,ax = plt.subplots()
# n_common = []
# n_unique = []
# from scipy.stats import mannwhitneyu,ks_2samp,ttest_ind,ttest_rel
# all_patients.pop(all_patients.index('Mel-16'))
# for p in all_patients:
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat/{}/combined/txt'.format(p))
#     pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     n_common_single = np.count_nonzero(pep['Proteins'].str.contains('common').values)
#     n_unique_single = np.count_nonzero(pep['Proteins'].str.contains('unique').values)
#     f_common = n_common_single / 114
#     f_unique = n_unique_single / 25188
#     n_common.append(f_common)
#     n_unique.append(f_unique)
# stats = ttest_ind(n_common,n_unique)
# print(stats)
# number = len(all_patients)
# ax.bar(x=[i for i in range(1,1+3*(number-1)+1,3)],height=n_common,label='common neoantigen')
# ax.bar(x=[i for i in range(2,2+3*(number-1)+1,3)],height=n_unique,label='unique neoantigen')
# xticks = []
# for ic,iu in zip([i for i in range(1,1+3*(number-1)+1,3)],[i for i in range(2,2+3*(number-1)+1,3)]):
#     xticks.append((ic+iu)/2)
# ax.set_xticks(xticks)
# ax.set_xticklabels(all_patients,fontsize=4,rotation=90)
# ax.set_xlabel('patients',fontsize=6)
# ax.set_ylabel('MS recovery rate',fontsize=6)
# ax.legend(loc='upper left',bbox_to_anchor=(1,1),frameon=False)
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# plt.savefig('MS_plot_concat.pdf',bbox_inches='tight')
# plt.close()

'''
plot concat stacked barplot
'''
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# fig,ax = plt.subplots()
# n_common = []
# n_unique = []
# from scipy.stats import mannwhitneyu,ks_2samp,ttest_ind,ttest_rel
# all_patients.pop(all_patients.index('Mel-16'))
# for p in all_patients:
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat/{}/combined/txt'.format(p))
#     pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     n_common_single = np.count_nonzero(pep['Proteins'].str.contains('common').values)
#     n_unique_single = np.count_nonzero(pep['Proteins'].str.contains('unique').values)
#     f_common = n_common_single / 114
#     f_unique = n_unique_single / 25188
#     n_common.append(f_common)
#     n_unique.append(f_unique)
# stats = ttest_ind(n_common,n_unique)
# print(stats)
# number = len(all_patients)
# ax.bar(x=np.arange(number),height=n_unique,label='unique neoantigen')
# ax.bar(x=np.arange(number),height=n_common,bottom=n_unique,label='common neoantigen')
# xticks = np.arange(number)
# ax.set_xticks(xticks)
# ax.set_xticklabels(all_patients,fontsize=4,rotation=90)
# ax.set_xlabel('patients',fontsize=6)
# ax.set_ylabel('MS recovery rate',fontsize=6)
# ax.legend(loc='upper left',bbox_to_anchor=(1,1),frameon=False)
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# plt.savefig('MS_plot_concat_stacked.pdf',bbox_inches='tight')
# plt.close()


'''
occurence, separate
'''
# dict_common = {}
# with open('result/MS_validation_common_unique2.fasta') as f:
#     for line in f:
#         if not line.startswith('>'):
#             dict_common[line.rstrip('\n')] = 0
# dict_unique = {}
# with open('result/MS_validation_unique_unique2.fasta') as f:
#     for line in f:
#         if not line.startswith('>'):
#             dict_unique[line.rstrip('\n')] = 0
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# all_patients.pop(all_patients.index('Mel-16'))
# for p in all_patients:
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common/{}/combined/txt'.format(p))
#     pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     for item in pep.index:
#         try:
#             dict_common[item] += 1
#         except KeyError:
#             continue
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_unique/{}/combined/txt'.format(p))
#     pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     for item in pep.index:
#         try:
#             dict_unique[item] += 1
#         except KeyError:
#             continue    
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# series_common = pd.Series(dict_common)
# series_common.name = 'occurence'
# series_common.to_csv('MS_common_occurence.txt',sep='\t')
# series_unique = pd.Series(dict_unique)
# series_unique.name = 'occurence'
# series_unique.to_csv('MS_unique_occurence.txt',sep='\t')
# fig,ax = plt.subplots()
# sns.ecdfplot(data=series_common.to_frame(),x='occurence',ax=ax)
# sns.ecdfplot(data=series_unique.to_frame(),x='occurence',ax=ax)
# import matplotlib.patches as mpatches
# ax.legend(handles=[mpatches.Patch(color=i) for i in ['#4B71B0','#DE8353']],labels=['common neoantigen','unique neoantigen'],loc='upper left',bbox_to_anchor=(1,1),frameon=False)
# #ax.text(x=0.5,y=0.1,s='Mann Whitney p-value: {}'.format(stats),transform=ax.transAxes)
# plt.savefig('MS_occurence.pdf',bbox_inches='tight');plt.close()


'''
occurence concat
'''
# dict_common = {}
# with open('result/MS_validation_common_unique2.fasta') as f:
#     for line in f:
#         if not line.startswith('>'):
#             dict_common[line.rstrip('\n')] = 0
# dict_unique = {}
# with open('result/MS_validation_unique_unique2.fasta') as f:
#     for line in f:
#         if not line.startswith('>'):
#             dict_unique[line.rstrip('\n')] = 0
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# all_patients.pop(all_patients.index('Mel-16'))
# for p in all_patients:
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat/{}/combined/txt'.format(p))
#     pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     for aa,protein in zip(*[pep.index,pep['Proteins']]):
#         if 'common' in protein:
#             try:
#                 dict_common[aa] += 1
#             except KeyError:
#                 continue
#         elif 'unique' in protein:
#             try:
#                 dict_unique[aa] += 1
#             except KeyError:
#                 continue  
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# series_common = pd.Series(dict_common)
# series_common.name = 'occurence'
# series_common.to_csv('MS_common_occurence.txt',sep='\t')
# series_unique = pd.Series(dict_unique)
# series_unique.name = 'occurence'
# series_unique.to_csv('MS_unique_occurence.txt',sep='\t')
# fig,ax = plt.subplots()
# sns.ecdfplot(data=series_common.to_frame(),x='occurence',ax=ax)
# sns.ecdfplot(data=series_unique.to_frame(),x='occurence',ax=ax)
# import matplotlib.patches as mpatches
# ax.legend(handles=[mpatches.Patch(color=i) for i in ['#4B71B0','#DE8353']],labels=['common neoantigen','unique neoantigen'],loc='upper left',bbox_to_anchor=(1,1),frameon=False)
# #ax.text(x=0.5,y=0.1,s='Mann Whitney p-value: {}'.format(stats),transform=ax.transAxes)
# plt.savefig('MS_occurence_concat.pdf',bbox_inches='tight');plt.close()

'''
occurence concat histogram or kde plot
'''
# dict_common = {}
# with open('result/MS_validation_common_unique2.fasta') as f:
#     for line in f:
#         if not line.startswith('>'):
#             dict_common[line.rstrip('\n')] = 0
# dict_unique = {}
# with open('result/MS_validation_unique_unique2.fasta') as f:
#     for line in f:
#         if not line.startswith('>'):
#             dict_unique[line.rstrip('\n')] = 0
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# all_patients.pop(all_patients.index('Mel-16'))
# for p in all_patients:
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat/{}/combined/txt'.format(p))
#     pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     for aa,protein in zip(*[pep.index,pep['Proteins']]):
#         if 'common' in protein:
#             try:
#                 dict_common[aa] += 1
#             except KeyError:
#                 continue
#         elif 'unique' in protein:
#             try:
#                 dict_unique[aa] += 1
#             except KeyError:
#                 continue  
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# series_common = pd.Series(dict_common)
# series_common.name = 'occurence'
# # series_common.to_csv('MS_common_occurence.txt',sep='\t')
# series_unique = pd.Series(dict_unique)
# series_unique.name = 'occurence'
# # series_unique.to_csv('MS_unique_occurence.txt',sep='\t')
# fig,ax = plt.subplots()
# sns.kdeplot(data=series_common.to_frame(),x='occurence',ax=ax,clip=[-0.05,15])
# sns.kdeplot(data=series_unique.to_frame(),x='occurence',ax=ax,clip=[-0.05,15])
# import matplotlib.patches as mpatches
# ax.legend(handles=[mpatches.Patch(color=i) for i in ['#4B71B0','#DE8353']],labels=['common neoantigen','unique neoantigen'],loc='upper left',bbox_to_anchor=(1,1),frameon=False)
# #ax.text(x=0.5,y=0.1,s='Mann Whitney p-value: {}'.format(stats),transform=ax.transAxes)
# plt.savefig('MS_occurence_concat_kde.pdf',bbox_inches='tight');plt.close()
# sys.exit('stop')


'''write peptide.txt to a merged xlsx file for publication, supp3 table'''
os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat')
all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
all_patients.pop(all_patients.index('Mel-16'))
df = pd.read_csv('result/frequency_stage3_verbosity1_uid_gene_coord_mean_mle.txt',sep='\t',index_col=0)
df['uid'] = [item.split(',')[1] for item in df.index]
df['aa'] = [item.split(',')[0] for item in df.index]
uid_2_ts_mean = pd.Series(data=df['tumor_specificity_mean'].values,index=df['uid'].values).to_dict()
uid_2_ts_mle = pd.Series(data=df['tumor_specificity_mle'].values,index=df['uid'].values).to_dict()
aa_2_uid = pd.Series(data=df['uid'].values,index=df['aa'].values).to_dict()
with pd.ExcelWriter('supp3_table.xlsx') as writer:
    for p in all_patients:
        os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat/{}/combined/txt'.format(p))
        pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
        col_mean = []
        col_mle = []
        for aa,protein in zip(*[pep.index,pep['Proteins']]):
            if pd.notna(protein):
                u = protein.split(',')[0].lstrip('>')
                mean = uid_2_ts_mean[u]
                mle = uid_2_ts_mle[u]
            else:
                mean,mle = pd.NA, pd.NA
            col_mean.append(mean)
            col_mle.append(mle)
        pep['tumor_specificity_mean'] = col_mean
        pep['tumor_specificity_mle'] = col_mle
        pep.to_excel(writer,sheet_name=p)
    os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
    # 18 neoantigens table
    df = pd.read_csv('MS_common_occurence.txt',sep='\t',index_col=0)
    df['uid'] = [aa_2_uid[item] for item in df.index]
    df['tumor_specificity_mean'] = [uid_2_ts_mean[item] for item in df['uid']]
    df['tumor_specificity_mle'] = [uid_2_ts_mle[item] for item in df['uid']]
    df.to_excel(writer,sheet_name='18_common_neoantigens')
    # add description
    dic = {
        'Sequence':'The amino acid sequence of the identified peptide',
        'N-term cleavage window':'Sequence window from -15 to 15 around the N-terminal cleavage site of this peptide',
        'C-term cleavage window':'Sequence window from -15 to 15 around the C-terminal cleavage site of this peptide',
        'Amino acid before':'The amino acid in the protein sequence before the peptide',
        'First amino acid':'The amino acid in the first position of the peptide sequence',
        'Second amino acid':'The amino acid in the second position of the peptide sequence',
        'Second last amino acid':'The amino acid in the second last position of the peptide sequence',
        'Last amino acid':'The amino acid in the last position of the peptide sequence',
        'Amino acid after':'The amino acid in the protein sequence after the peptide',
        'N Count':'The number of instances of the "N" amino acid contained within the sequence, N indicates amino acid letter',
        'Length':'The length of the sequence stored in the column "Sequence"',
        'Missed cleavages':'Number of missed enzymatic cleavages',
        'Mass':'Monoisotopic mass of the peptide',
        'Proteins':'Identifiers of proteins this peptide is associated with',
        'Leading razor protein':'Identifier of the leading protein in the protein group which uses this peptide for quantification. (Either unique or razor)',
        'Start position':'Position of the first amino acid of this peptide in the protein sequence. (one-based)',
        'End position':'Position of the last amino acid of this peptide in the protein sequence. (one-based)',
        'Unique (Groups)':'When marked with "+", this particular peptide is unique to a single protein group in the proteinGroups file',
        'Unique (Proteins)':'When marked with "+", this particular peptide is unique to a single protein sequence in the fasta file(s)',
        'Charges':'All charge states that have been observed',
        'PEP':'Posterior Error Probability of the identification. This value essentially operates as a p-value, where smaller is more significant',
        'Score':'Highest Andromeda score for the associated MS/MS spectra',
        'Experient [n]': 'Number of evidence entries for this "Experiment [n]"',
        'Intensity':'Summed up eXtracted Ion Current (XIC) of all isotopic clusters associated with the identified AA sequence. In case of a labeled experiment this is the total intensity of all the isotopic patterns in the label cluster',
        'Reverse':'When marked with "+", this particular peptide was found to be part of a protein derived from the reversed part of the decoy database. These should be removed for further data analysis',
        'Potential contaminant':'When marked with '+', this particular peptide was found to be part of a commonly occurring contaminant. These should be removed for further data analysis',
        'id':'A unique (consecutive) identifier for each row in the peptides table, which is used to cross-link the information in this table with the information stored in the other tables',
        'Protein group IDs':'The identifiers of the protein groups this peptide was linked to, referenced against the proteinGroups table',
        'Mod. peptide IDs':'Identifier(s) for peptide sequence(s), associated with the peptide, referenced against the corresponding modified peptides table',
        'Evidence IDs':'Identifier(s) for analyzed peptide evidence associated with the protein group referenced against the evidence table',
        'MS/MS IDs':'The identifiers of the MS/MS scans identifying this peptide, referenced against the msms table',
        'Best MS/MS':'The identifier of the best (in terms of quality) MS/MS scan identifying this peptide, referenced against the msms table',
        'Oxidation (M) site IDs':'Identifier(s) for site(s) associated with the protein group, which show(s) evidence of the modification, referenced against the appropriate modification site file',
        'Taxonomy IDs':'Taxonomy identifiers',
        'MS/MS Count':'The number of MS/MS evidence',
        'tumor_specificity_mean':'The average raw read counts for each parental NeoJunction',
        'tumor_specificity_mle':'The Tumor specificity score derived using Maximum Likelihood Estimation(MLE)',
    }
    pd.Series(data=dic,name='description').to_excel(writer,sheet_name='description')


sys.exit('stop')





## tumorigenesis
# df = pd.read_csv('result/toppfun_stage3_0.6.txt',sep='\t',encoding = 'unicode_escape')
# df = df.loc[df['Category']=='GO: Biological Process',:]
# df['convert'] = -np.log10(df['q-value FDR B&H'].values)
# df = df.iloc[:5,:]
# df.sort_values(by='convert',axis=0,ascending=True,inplace=True)
# fig,ax = plt.subplots()
# ax.barh(y=np.arange(df.shape[0]),width=[item for item in df['convert']],color='#FF9A91')
# ax.set_yticks(np.arange(df.shape[0]))
# ax.set_yticklabels([item for item in df['Name']])
# ax.set_title('GO: Biological Process')
# ax.set_xlabel('-Log10(adjusted_pval)')
# plt.savefig('go_biological_process.pdf',bbox_inches='tight')
# plt.close()

# df = pd.read_csv('result/toppfun_stage3_0.6.txt',sep='\t',encoding = 'unicode_escape')
# df = df.loc[df['Category']=='GO: Cellular Component',:]
# df['convert'] = -np.log10(df['q-value FDR B&H'].values)
# df = df.iloc[:5,:]
# df.sort_values(by='convert',axis=0,ascending=True,inplace=True)
# fig,ax = plt.subplots()
# ax.barh(y=np.arange(df.shape[0]),width=[item for item in df['convert']],color='#FF9A91')
# ax.set_yticks(np.arange(df.shape[0]))
# ax.set_yticklabels([item for item in df['Name']])
# ax.set_title('GO: Cellular Component')
# ax.set_xlabel('-Log10(adjusted_pval)')
# plt.savefig('go_cellular_component.pdf',bbox_inches='tight')
# plt.close()

# df = pd.read_csv('result/toppfun_stage3_0.6.txt',sep='\t',encoding = 'unicode_escape')
# df = df.loc[df['Category']=='Pathway',:]
# df['convert'] = -np.log10(df['q-value FDR B&H'].values)
# df = df.iloc[:5,:]
# df.sort_values(by='convert',axis=0,ascending=True,inplace=True)
# fig,ax = plt.subplots()
# ax.barh(y=np.arange(df.shape[0]),width=[item for item in df['convert']],color='#FF9A91')
# ax.set_yticks(np.arange(df.shape[0]))
# ax.set_yticklabels([item for item in df['Name']])
# ax.set_title('Pathway')
# ax.set_xlabel('-Log10(adjusted_pval)')
# plt.savefig('pathway.pdf',bbox_inches='tight')
# plt.close()



### Step2: necessary secondary results
# df = pd.read_csv('./result/frequency_stage3_verbosity1_uid.txt',sep='\t',index_col=0)
# snaf.downstream.add_gene_symbol_frequency_table(df=df).to_csv('./result/frequency_stage3_verbosity1_uid_gene_symbol.txt',sep='\t')
# jcmq = snaf.JunctionCountMatrixQuery.deserialize(name='./result/after_prediction.p')
# jcmq.visualize(uid='ENSG00000167291:E38.6-E39.1',sample='TCGA-DA-A1I1-06A-12R-A18U-07.bed',outdir='./result')

### Step3: downstream analysis (patient level and neoantigen level)

'''patient analysis'''
# 1. survival analysis
# survival = pd.read_csv('TCGA-SKCM.survival.tsv',sep='\t',index_col=0)  # 463
# burden = pd.read_csv('result/burden_stage0.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
# burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
# burden_output,quantiles = snaf.survival_analysis(burden,survival,n=2,stratification_plot='result/stage0_stratify.pdf',survival_plot='result/stage0_survival.pdf')
# burden_output.to_csv('result/to_nathan_stage0_neojunction_encode.txt',sep='\t')



# 2. mutation analysis
# mutation = pd.read_csv('TCGA-SKCM.mutect2_snv.tsv',sep='\t',index_col=0)  # 467 samples have mutations
# mutation = mutation.loc[mutation['filter']=='PASS',:]
# burden = pd.read_csv('result/burden_stage3.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
# burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
# snaf.mutation_analysis(mode='compute',burden=burden,mutation=mutation,output='result/stage3_mutation.txt')


'''neoantigen analysis'''
# 1. physicalchemical properties and occurence frequency
# snaf.downstream.analyze_neoantigens(freq_path='result/frequency_stage3_verbosity1_uid.txt',junction_path='result/burden_stage0.txt',total_samples=472,outdir='result',mers=None,fasta=False)
# snaf.run_dash_T_antigen(input_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/result/shared_vs_unique_neoantigen_all.txt',
#                         output_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/result')
# snaf.downstream.plot_umap_neoantigen(df_path='result/mer9_umap_embed_df.txt',outdir='result')
        




























