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

# 4 common neoantigens for immuno assay
# cand = pd.read_csv('result_new/T_candidates/T_antigen_candidates_all.txt',sep='\t',index_col=0)
# outdir = '4_common_inspect'
# uid = 'ENSG00000071991:E12.1-E13.2_66509195'
# criterion=[('netMHCpan_el', 0, '<=', 2),('deepimmuno_immunogenicity',1,'==','True')]
# snaf.gtex_visual_combine_plotly(uid=uid,norm=False,outdir=outdir,tumor=df)
# snaf.gtex_visual_combine_plotly(uid=uid,norm=True,outdir=outdir,tumor=df)
# snaf.JunctionCountMatrixQuery.deserialize('result_new/after_prediction.p').visualize(uid,'TCGA-FS-A4FB-06A-11R-A266-07.bed',outdir=outdir,tumor=False,criterion=criterion)
# sub_cand = cand.loc[cand['uid']==uid,:]
# sub_cand.to_csv('{}/{}_cand.txt'.format(outdir,uid.replace(':','_')),sep='\t')

# whether there's a HLA-C bias in prediction
# cand = pd.read_csv('result_new/T_candidates/T_antigen_candidates_all.txt',sep='\t',index_col=0)
# uid = [p+','+h for p,h in zip(cand['peptide'],cand['hla'])]
# cand['key'] = uid 
# cand = cand.loc[cand['key'].duplicated(),:]
# col = 'immunogenicity'
# hla_c_bind = cand.loc[cand['hla'].str.startswith('HLA-C'),col]
# hla_b_bind = cand.loc[cand['hla'].str.startswith('HLA-B'),col]
# hla_a_bind = cand.loc[cand['hla'].str.startswith('HLA-A'),col]
# fig,ax = plt.subplots()
# sns.barplot(data=[[[len(hla_a_bind)]],[len(hla_b_bind)],[len(hla_c_bind)]])
# ax.set_xticklabels(['hla-a','hla-b','hla-c'])
# ax.set_ylabel('number of hits for each major allele')
# plt.savefig('c_bias_num.pdf',bbox_inches='tight')
# plt.close()
# fig,ax = plt.subplots()
# sns.violinplot(data=[hla_a_bind,hla_b_bind,hla_c_bind],ax=ax)
# ax.set_xticklabels(['hla-a','hla-b','hla-c'])
# ax.set_ylabel(col)
# plt.savefig('c_bias_value.pdf',bbox_inches='tight')
# plt.close()

'''T cell neoantigen'''
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=30,add_control=add_control,outdir='result_new',filter_mode='maxmin')
# sample_to_hla = pd.read_csv('sample_hla.txt',sep='\t',index_col=0)['hla'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq.run(hlas=hlas,outdir='./result_new')
# snaf.JunctionCountMatrixQuery.generate_results(path='./result_new/after_prediction.p',outdir='./result_new')

# do survival and mutation analysis
# survival = pd.read_csv('TCGA-SKCM.survival.tsv',sep='\t',index_col=0)
# for stage in [0,2,3]:
#     burden_df = pd.read_csv('result_new/burden_stage{}.txt'.format(stage),sep='\t',index_col=0)
#     burden_df.rename(columns=lambda x:'-'.join(x.split('-')[:4]),inplace=True)
#     burden_output,_ = snaf.survival_analysis(burden_df,survival,2,stratification_plot='result_new/survival/stage{}_stratify.pdf'.format(stage),
#                                              survival_plot='result_new/survival/stage{}_survival.pdf'.format(stage))
#     if stage == 3:
#         burden_output.to_csv('result_new/survival/burden3_patient_high_low_group.txt',sep='\t')

# mutation = pd.read_csv('TCGA-SKCM.mutect2_snv.tsv',sep='\t',index_col=0)
# burden3 = pd.read_csv('result_new/burden_stage3.txt',sep='\t',index_col=0)
# burden3.rename(columns=lambda x:'-'.join(x.split('-')[:4]),inplace=True)
# snaf.mutation_analysis(mode='compute',burden=burden3,mutation=mutation,output='result_new/survival/mutation.txt')
# snaf.mutation_analysis(mode='plot',burden=burden3,mutation=mutation,output='result_new/survival/CAMKK2_mutation.txt',genes_to_plot=['CAMKK2'])

# snaf.downstream.survival_regression(freq='result_new/frequency_stage3_verbosity1_uid_gene_symbol_coord_mean_mle.txt',remove_quote=True,
#                                     rename_func=lambda x:'-'.join(x.split('-')[:4]),survival='TCGA-SKCM.survival.tsv',
#                                     pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',outdir='result_new/survival',mode='binary')

# do DEG and GO
# snaf.downstream.prepare_DEG_analysis('result_new/burden_stage3.txt','result_new/survival/burden3_patient_high_low_group.txt',
#                                      rename_func=lambda x:'-'.join(x.split('-')[:4]),
#                                      outdir='result_new/survival',encoding={'low':'1','high':'2'})
# snaf.downstream.visualize_DEG_result('result_new/survival/DEGs-LogFold_0.0_adjp/GE.low_vs_high.txt',up_cutoff=0.58,down_cutoff=-0.58,
#                                      mode='static',outdir='result_new/survival',genes_to_highlight=['LST1','HCST','IL32','CD3D','S100A8','MZB1','IGLC4','ADAM10','ARFGEF2','MIB1','KIF3B','TNPO1','PTPN11','ANKRD52','TGFBR1'])

# visualize RBP specifically
# rbps = ['YTHDC2', 'SDAD1', 'TARDBP', 'UCHL5', 'SRSF1', 'PUM1', 'GNL3', 'PNPT1', 'DHX30', 'GRSF1', 'HLTF', 'PSIP1', 'DDX52', 'EIF4B', 'SND1', 'PABPC1', 'WDR43', 'AGO2', 'STAU1', 'SSB', 'HNRNPU', 'EIF2S2', 'FAM120A', 'RBM39', 'SNRNP200', 'EIF4G1', 'HSPD1', 'DDX55', 'DDX21', 'PPIG', 'DKC1', 'FIP1L1', 'LARP4', 'YWHAG', 'RCC2', 'EIF4G2', 'MSI2', 'WDR3', 'EIF3A', 'PRPF4', 'PRPF8', 'SERBP1', 'CELF1', 'AQR', 'SUPV3L1', 'FTO', 'NUSAP1', 'XPO5', 'AGO1', 'FMR1', 'SF3B1', 'DNAJC21', 'TROVE2', 'DDX27', 'ADAR', 'KHDRBS1', 'TIAL1', 'CPEB4', 'KRR1', 'FXR1', 'AGO3', 'XPO1', 'NIPBL', 'EIF2S1', 'HNRNPF', 'XRN1', 'DDX42', 'ZC3H11A', 'PUM2', 'AKAP1', 'DROSHA', 'EFTUD2', 'MTPAP', 'DDX1', 'NUFIP2', 'IGF2BP2', 'G3BP1', 'QKI', 'IGF2BP3', 'LARP7', 'GEMIN5', 'DDX3X', 'DDX6', 'TIA1', 'SLTM', 'TRIM56', 'MBNL1', 'SRFBP1', 'AGGF1', 'SUGP2', 'FUBP1', 'ZNF800', 'SF3A3', 'G3BP2', 'MAK16', 'RBM25', 'METAP2', 'CPSF6', 'ZRANB2', 'SSRP1', 'NOLC1', 'UTP18', 'NSUN2', 'FASTKD2', 'WRN', 'ESF1', 'ILF3', 'RBFOX2', 'SART3', 'PABPC4', 'PHF6', 'ABCF1', 'NAA15', 'SRP68', 'XRN2', 'MATR3']
# snaf.downstream.visualize_DEG_result('result_new/survival/DEGs-LogFold_0.0_adjp/GE.low_vs_high.txt',up_cutoff=0.58,down_cutoff=-0.58,
#                                      mode='static',outdir='result_new/survival',genes_to_highlight=rbps)

# snaf.downstream.prepare_GO_analysis('result_new/survival/DEGs-LogFold_0.0_adjp/GE.low_vs_high.txt',outdir='result_new/survival',direction='<',lc_cutoff=-0.58,adjp_cutoff=0.05)
# snaf.downstream.visualize_GO_result(path_list=['result_new/survival/GO_Elite_result_GeneOntology/GO-Elite_results/CompleteResults/ORA/archived-20221201-205821/gene_list-GO.txt','result_new/survival/GO_Elite_result_BioMarkers/GO-Elite_results/CompleteResults/ORA/archived-20221201-205720/gene_list-BioMarkers.txt'],
#                                     skiprows_list=[17,16],category_list=['Ontology Name','Gene-Set Name'],outdir='result_new/survival',
#                                     mode='static',ontology_to_highlight={'Adult Peripheral Blood Activated T cell (PMID32214235 top 100)':'T cells','antigen binding':'antigen binding','complement activation':'Complement Activation','immune response':'immune response','humoral immune response':'humoral immune response'},ylims=(10e-50,10e-1))
# snaf.downstream.visualize_GO_result(path_list=['result_new/survival/GO_Elite_result_GeneOntology/GO-Elite_results/CompleteResults/ORA/archived-20230407-100025/gene_list-GO.txt'],
#                                     skiprows_list=[17],category_list=['Ontology Name'],outdir='result_new/survival',
#                                     mode='static',ontology_to_highlight={'metabolic process':'metabolic process','nucleic acid metabolic process':'nucleic acid metabolic process','nitrogen compound metabolic process':'nitrogen compound metabolic process','regulation of cellular metabolic process':'regulation of cellular metabolic process','response to DNA damage stimulus':'response to DNA damage stimulus'},ylims=(10e-80,10e-1))

# # coverage analysis
# snaf.downstream.get_coverage(t_result='result_new/T_candidates/T_antigen_candidates_all.txt',allele='A')
# sys.exit('stop')


'''B cell neoantigen'''
# # long-read
# membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df,add_control=add_control,not_in_db=False,outdir='result_new/surface',filter_mode='maxmin')
# surface.run(uids=membrane_tuples,outdir='result_new/surface',prediction_mode='long_read',n_stride=2,
#             gtf='/data/salomonis2/LabFiles/Frank-Li/refactor/data/2021UHRRIsoSeq_SQANTI3_filtered.gtf',
#             tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',serialize=True)
# surface.generate_full_results(outdir='result_new/surface',freq_path='result_new/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt',mode='long_read',validation_gtf=None)

# # short-read
# membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df,add_control=add_control,not_in_db=False,outdir='result_new/surface',filter_mode='maxmin')
# surface.run(uids=membrane_tuples,outdir='result_new/surface',prediction_mode='short_read',n_stride=2,
#             gtf=None,
#             tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm',serialize=True)
# surface.generate_full_results(outdir='result_new/surface',freq_path='result_new/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt',mode='short_read',validation_gtf='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/SQANTI-all/collapse_isoforms_classification.filtered_lite.gtf')

# surface.run_dash_B_antigen(pkl='result_new/surface/surface_antigen_lr.p',candidates='result_new/surface/candidates_3_lr_None_False.txt',prediction_mode='long_read',
#                            python_executable='/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7')




## MS validation
# df = pd.read_csv('result_new/frequency_stage3_verbosity1_uid_gene_symbol_coord_mean_mle.txt',sep='\t',index_col=0)
# df_common = df.loc[df['n_sample']>71,:]
# df_unique = df.loc[df['n_sample']==1,:]
# with open('result_new/MS_validation_common.fasta','w') as f:
#     for row in df_common.itertuples():
#         peptide, uid = row.Index.split(',')
#         f.write('>{}\n{}\n'.format(uid,peptide))
# with open('result_new/MS_validation_unique.fasta','w') as f:
#     for row in df_unique.itertuples():
#         peptide, uid = row.Index.split(',')
#         f.write('>{}\n{}\n'.format(uid,peptide))
# for f in ['MS_validation_common','MS_validation_unique']:
#     snaf.proteomics.remove_redundant('./result_new/{}.fasta'.format(f),'./result_new/{}_unique.fasta'.format(f))
#     snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot_9_10_mers_unique.fasta',
#                                       fa2_path='./result_new/{}_unique.fasta'.format(f),outdir='./result_new',write_unique2=True,prefix='{}_'.format(f))


# with open('result_new/MS_validation_common_unique2.fasta','r') as f1, open('result_new/MS_validation_unique_unique2.fasta','r') as f2, open('result_new/MS_validation_concat_unique2.fasta','w') as f3:
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

# # do some move
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat')
# all_patients = subprocess.run('for folder in *; do echo $folder; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# for p in all_patients:
#     print(p)
#     os.mkdir('/data/salomonis-archive/MS/melanoma/raw/{}'.format(p))
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat/{}'.format(p))
#     des_folder = '/data/salomonis-archive/MS/melanoma/raw/{}'.format(p)
#     subprocess.run('for file in *.raw; do mv $file {}; done'.format(des_folder),shell=True)
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat')


# db_concat = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/MS_validation_concat_unique2.fasta']
# os.chdir('/data/salomonis-archive/MS/melanoma/raw')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# for p in all_patients:
#     os.chdir('/data/salomonis-archive/MS/melanoma/raw/{}'.format(p))
#     all_samples = subprocess.run('for file in *.raw; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     inputs = []
#     pwd = os.getcwd()
#     for s in all_samples:
#         inputs.append(os.path.join(pwd,s))
#     outdir='/data/salomonis-archive/MS/melanoma/raw/{}'.format(p)
#     snaf.proteomics.set_maxquant_configuration(dbs=db_concat,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,outdir=outdir)


'''
plot concat stacked barplot
'''
# os.chdir('/data/salomonis-archive/MS/melanoma/raw')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
# fig,ax = plt.subplots()
# n_common = []
# n_unique = []
# from scipy.stats import mannwhitneyu,ks_2samp,ttest_ind,ttest_rel
# all_patients.pop(all_patients.index('Mel-16'))
# for p in all_patients:
#     os.chdir('/data/salomonis-archive/MS/melanoma/raw/{}/combined/txt'.format(p))
#     pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     n_common_single = np.count_nonzero(pep['Proteins'].str.contains('common').values)
#     n_unique_single = np.count_nonzero(pep['Proteins'].str.contains('unique').values)
#     f_common = n_common_single / 613
#     f_unique = n_unique_single / 16753
#     n_common.append(f_common)
#     n_unique.append(f_unique)
# stats = ttest_rel(n_common,n_unique)
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
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
# plt.savefig('MS_plot_concat_stacked.pdf',bbox_inches='tight')
# plt.close()


'''
plot concat side-by-side barplot
'''
# os.chdir('/data/salomonis-archive/MS/melanoma/raw')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
# fig,ax = plt.subplots()
# n_common = []
# n_unique = []
# from scipy.stats import mannwhitneyu,ks_2samp,ttest_ind,ttest_rel
# all_patients.pop(all_patients.index('Mel-16'))
# for p in all_patients:
#     os.chdir('/data/salomonis-archive/MS/melanoma/raw/{}/combined/txt'.format(p))
#     pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     n_common_single = np.count_nonzero(pep['Proteins'].str.contains('common').values)
#     n_unique_single = np.count_nonzero(pep['Proteins'].str.contains('unique').values)
#     f_common = n_common_single / 613
#     f_unique = n_unique_single / 16753
#     n_common.append(f_common)
#     n_unique.append(f_unique)
# stats = mannwhitneyu(n_common,n_unique)   # should use related as they are related, p=0.006
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
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
# plt.savefig('MS_plot_concat.pdf',bbox_inches='tight')
# plt.close()

'''
occurence concat
'''
# dict_common = {}
# with open('result_new/MS_validation_common_unique2.fasta') as f:
#     for line in f:
#         if not line.startswith('>'):
#             dict_common[line.rstrip('\n')] = 0
# dict_unique = {}
# with open('result_new/MS_validation_unique_unique2.fasta') as f:
#     for line in f:
#         if not line.startswith('>'):
#             dict_unique[line.rstrip('\n')] = 0
# os.chdir('/data/salomonis-archive/MS/melanoma/raw')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
# all_patients.pop(all_patients.index('Mel-16'))
# for p in all_patients:
#     os.chdir('/data/salomonis-archive/MS/melanoma/raw/{}/combined/txt'.format(p))
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
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
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
# plt.savefig('MS_occurence_concat.pdf',bbox_inches='tight');plt.close()


'''
occurence concat histogram or kde plot
'''
# dict_common = {}
# with open('result_new/MS_validation_common_unique2.fasta') as f:
#     for line in f:
#         if not line.startswith('>'):
#             dict_common[line.rstrip('\n')] = 0
# dict_unique = {}
# with open('result_new/MS_validation_unique_unique2.fasta') as f:
#     for line in f:
#         if not line.startswith('>'):
#             dict_unique[line.rstrip('\n')] = 0
# os.chdir('/data/salomonis-archive/MS/melanoma/raw')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
# all_patients.pop(all_patients.index('Mel-16'))
# for p in all_patients:
#     os.chdir('/data/salomonis-archive/MS/melanoma/raw/{}/combined/txt'.format(p))
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
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
# series_common = pd.Series(dict_common)
# series_common.name = 'occurence'
# series_unique = pd.Series(dict_unique)
# series_unique.name = 'occurence'
# fig,ax = plt.subplots()
# sns.kdeplot(data=series_common.to_frame(),x='occurence',ax=ax,clip=[-0.05,15])
# sns.kdeplot(data=series_unique.to_frame(),x='occurence',ax=ax,clip=[-0.05,15])
# import matplotlib.patches as mpatches
# ax.legend(handles=[mpatches.Patch(color=i) for i in ['#4B71B0','#DE8353']],labels=['common neoantigen','unique neoantigen'],loc='upper left',bbox_to_anchor=(1,1),frameon=False)
# plt.savefig('MS_occurence_concat_kde.pdf',bbox_inches='tight');plt.close()

'''
further explore the identity of common neoantigens
'''
# identified_common_neoantigens = []
# os.chdir('/data/salomonis-archive/MS/melanoma/raw')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
# all_patients.pop(all_patients.index('Mel-16'))
# for p in all_patients:
#     os.chdir('/data/salomonis-archive/MS/melanoma/raw/{}/combined/txt'.format(p))
#     pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#     pep = pep.loc[pep['Proteins'].notna(),:]
#     for aa,protein,pval,s in zip(*[pep.index,pep['Proteins'],pep['PEP'],pep['Score']]):
#         if 'common' in protein:
#             data = (aa,protein.lstrip('>').split(',')[0],pval,s,p)
#             identified_common_neoantigens.append(data)
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
# pd.DataFrame.from_records(data=identified_common_neoantigens,columns=['peptide','uid','PEP','Score','patient']).to_csv('identified_common_neoantigens.txt',sep='\t')

'''
inspect those candidates for junction validity and tumor specificity
'''
# control_bam_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A208-11A-51R-A157-07.bam',
#                  '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-Kidney/KIRC/TCGA-CJ-5680-11A-01R-1541-07.bam',
#                  '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-STAD/TCGA-BR-6453-11A-01R-1802-13.bam',
#                  '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-ESCA/TCGA-L5-A4OG-11A-12R-A260-31.bam',
#                  '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-THCA/TCGA-EL-A3ZP-11A.bam',
#                  '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A0H7-11A-13R-A089-07.bam']
# control_bai_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A208-11A-51R-A157-07.bam.bai',
#                  '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-Kidney/KIRC/TCGA-CJ-5680-11A-01R-1541-07.bam.bai',
#                  '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-STAD/TCGA-BR-6453-11A-01R-1802-13.bam.bai',
#                  '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-ESCA/TCGA-L5-A4OG-11A-12R-A260-31.bam.bai',
#                  '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-THCA/TCGA-EL-A3ZP-11A.bam.bai',
#                  '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A0H7-11A-13R-A089-07.bam.bai']

# control_bam_path_list = ['/data/salomonis-archive/BAMs/PublicDatasets/E-MTAB-2836-Grch38_Deep-Healthy-PanTissue/ERR315409_1.bam']
# control_bai_path_list = ['/data/salomonis-archive/BAMs/PublicDatasets/E-MTAB-2836-Grch38_Deep-Healthy-PanTissue/ERR315409_1.bam.bai']

control_bam_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A208-11A-51R-A157-07.bam']
control_bai_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A208-11A-51R-A157-07.bam.bai']

def unit_run(uid,base_sample,region,task_name):
    uid = uid
    sample = base_sample + '.bed'
    region = region
    bam_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-SKCM/TCGA_SKCM-BAMs/bams/{}.bam'.format(base_sample)] + control_bam_path_list
    bai_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-SKCM/TCGA_SKCM-BAMs/bams/{}.bam.bai'.format(base_sample)] + control_bai_path_list
    sif_anno_path = '/data/salomonis2/software/ggsashimi'
    outdir = 'Frank_inspection/sashimi'
    bam_contig_rename = [False,False]
    criterion=[('netMHCpan_el', 0, '<=', 2),('deepimmuno_immunogenicity',1,'==','True')]

    # snaf.gtex_visual_combine_plotly(uid=uid,outdir='Frank_inspection',norm=False,tumor=df)
    # snaf.gtex_visual_combine_plotly(uid=uid,outdir='Frank_inspection',norm=True,tumor=df)
    # snaf.JunctionCountMatrixQuery.deserialize('result_new/after_prediction.p').visualize(uid=uid,sample=sample,outdir='Frank_inspection',criterion=criterion)
    dff = snaf.gtex_visual_combine(uid=uid,outdir='Frank_inspection',norm=False,tumor=df)
    dff.to_csv('HAAASFETL_adata_df.txt',sep='\t')
    # snaf.prepare_sashimi_plot(bam_path_list,bai_path_list,outdir,sif_anno_path,bam_contig_rename,query_region=region,skip_copy=False, min_junction=1,task_name=task_name)

def flank_chrom(chrom,offset):
    chrom = chrom.split('(')[0]
    sec1 = chrom.split(':')[0]
    sec2 = chrom.split(':')[1].split('-')[0]
    sec3 = chrom.split(':')[1].split('-')[1]
    new_sec2 = int(sec2) - offset[0]
    new_sec3 = int(sec3) + offset[1]
    assemble = '{}:{}-{}'.format(sec1,new_sec2,new_sec3)
    return assemble



# df = pd.read_csv('candidates.csv',sep=',')
# offset = (1000,1000)
# lis = [(row.UID,row.sample,flank_chrom(row.chrom,offset),row.Index) for row in df.itertuples()]

# item = lis[48]
# print(item[0])
# unit_run(item[0],item[1],item[2],'{}_offset_{}'.format(item[3],offset))

# specifically, generate publication-quality sashimi for 8 validated hits
'''
VAPGEAKNL
ENSG00000185989:E5.1_114128541-E18.1
Mel4
TCGA-DA-A1HY-06A-11R-A18T-07
RASA3
gene is highly expressed, but splicing is unique, very tumor specific
adipose tissue ERR315343 use this as control instead breast is also highly expressed TCGA-BH-A0H7-11A-13R-A089-07
'''
# unit_run(uid='ENSG00000185989:E5.1_114128541-E18.1',
#          base_sample='TCGA-DA-A1HY-06A-11R-A18T-07',
#          region=flank_chrom('chr13:114073837-114128541',(10000,10000)),
#          task_name='_VAPGEAKNL')

'''
YALANIKWI
ENSG00000152760:E7.2-I7.1
Mel12
TCGA-FS-A1ZF-06A-12R-A18S-07
DYNLT5
fallopian tube,ERR579134
'''
# unit_run(uid='ENSG00000152760:E7.2-I7.1',
#          base_sample='TCGA-FS-A1ZF-06A-12R-A18S-07',
#          region=flank_chrom('chr1:66775578-66775579',(1000,5000)),
#          task_name='_YALANIKWI')

'''
HAAASFETL
ENSG00000132185:E2.3-E2.5
Mel15
TCGA-EB-A3HV-01A-11R-A21D-07
FCRLA
tonsil, ERR579133, super tumor specific
end up using TCGA-BR-6453-11A-01R-1802-13 stomatch instead
'''
unit_run(uid='ENSG00000132185:E2.3-E2.5',
         base_sample='TCGA-EB-A3HV-01A-11R-A21D-07',
         region=flank_chrom('chr1:161710492-161710760',(5000,1000)),
         task_name='_HAAASFETL')
sys.exit('stop')

'''
TELQRTLSL
ENSG00000151092:E23.1-E24.1
Mel20/26
TCGA-EB-A44O-01A-11R-A266-07
NGLY1
end up using TCGA-BH-A0H7-11A-13R-A089-07
'''
# unit_run(uid='ENSG00000151092:E23.1-E24.1',
#          base_sample='TCGA-EB-A44O-01A-11R-A266-07',
#          region=flank_chrom('chr3:25736400-25737334',(3000,1000)),
#          task_name='_TELQRTLSL')

'''
KEKLDQLVY
ENSG00000100225:E8.2-E9.2_32491143
Mel27
TCGA-ER-A2NE-06A-21R-A18T-07
FBXO7
end up using TCGA-BH-A208-11A-51R-A157-07
'''
unit_run(uid='ENSG00000100225:E8.2-E9.2_32491143',
         base_sample='TCGA-ER-A2NE-06A-21R-A18T-07',
         region=flank_chrom('chr22:32487828-32491143',(500,500)),
         task_name='_KEKLDQLVY')
sys.exit('stop')


'''compare against the proliferative tissues'''
# # get reduced junction
# df = snaf.get_reduced_junction_matrix(pc='counts.TCGA-SKCM.txt',pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')

# # run SNAF but just get the maxmin text file
# netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
# db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
# tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
# gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
# proliferate_db = snaf.remove_trailing_coord('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/proliferating/counts.original.txt')
# proliferate_db.rename(columns=lambda x:x+'_proliferate',inplace=True)
# add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db,'proliferate':proliferate_db}

# snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
# surface.initialize(db_dir=db_dir)
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=30,add_control=add_control,outdir='result_proliferate',filter_mode='maxmin')

'''
original shape: (951494, 472)
current shape: (101974, 472)
2023-04-07 07:48:40 starting initialization
Current loaded gtex cohort with shape (100039, 2629)
Adding cohort tcga_control with shape (101974, 705) to the database
now the shape of control db is (101974, 3334)
Adding cohort gtex_skin with shape (90788, 313) to the database
now the shape of control db is (101974, 3647)
Adding cohort proliferate with shape (88752, 169) to the database
now the shape of control db is (101974, 3816)
2023-04-07 07:50:29 finishing initialization
2023-04-07 07:50:29 starting surface antigen initialization
2023-04-07 07:50:53 finished surface antigen initialization
reduce valid NeoJunction from 101974 to 22936 because they are present in GTEx
reduce valid Neojunction from 22936 to 17670 because they are present in added control tcga_control
reduce valid Neojunction from 17670 to 16799 because they are present in added control gtex_skin
reduce valid Neojunction from 16799 to 15026 because they are present in added control proliferate
'''

# get gene_list that are removed in last step
# filter_df = pd.read_csv('result_proliferate/NeoJunction_statistics_maxmin.txt',sep='\t',index_col=0)
# filter_df_removed = filter_df.loc[(filter_df['cond']==True)&(filter_df['cond_add_tcga_control']==True)&(filter_df['cond_add_gtex_skin']==True)&(filter_df['cond_add_proliferate']==False),:]
# filter_df_removed.to_csv('result_proliferate/filter_df_filtered.txt',sep='\t')
# df = filter_df_removed.loc[:,['max','mean_add_proliferate']]
# sns.stripplot(data=df,orient='h')
# plt.savefig('result_proliferate/boxplot.pdf',bbox_inches='tight')
# plt.close()
# genes = list(set([item.split(':')[0] for item in filter_df_removed.index]))
# symbols = snaf.downstream.ensemblgene_to_symbol(query=genes,species='human')
# with open('result_proliferate/gene_list.txt','w') as f:
#     for s in symbols:
#         if s != 'unknown_gene':
#             f.write('{}\n'.format(s))
# snaf.visualize_GO_result(path_list=['/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_proliferate/GO_Elite_result_GeneOntology/GO-Elite_results/CompleteResults/ORA/archived-20230407-085012/gene_list-GO.txt'],
#                          skiprows_list=[17],category_list=['Ontology Name'],outdir='result_proliferate',mode='static',ylims=(10e-25,10e-1),
#                          ontology_to_highlight=
#                              {'cell cycle':'cell cycle',
#                               'cell division':'cell division',
#                               'chromosome segregation':'chromosome segregation',
#                               'spindle':'spindle',
#                               'cytoskeleton':'cytoskeleton',
#                               'regulation of cell cycle process':'regulation of cell cycle process'})


'''write peptide.txt to a merged xlsx file for publication, supp3 table'''
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
# all_patients.pop(all_patients.index('Mel-16'))
# df = pd.read_csv('result/frequency_stage3_verbosity1_uid_gene_coord_mean_mle.txt',sep='\t',index_col=0)
# df['uid'] = [item.split(',')[1] for item in df.index]
# df['aa'] = [item.split(',')[0] for item in df.index]
# uid_2_ts_mean = pd.Series(data=df['tumor_specificity_mean'].values,index=df['uid'].values).to_dict()
# uid_2_ts_mle = pd.Series(data=df['tumor_specificity_mle'].values,index=df['uid'].values).to_dict()
# aa_2_uid = pd.Series(data=df['uid'].values,index=df['aa'].values).to_dict()
# with pd.ExcelWriter('supp3_table.xlsx') as writer:
#     for p in all_patients:
#         os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_concat/{}/combined/txt'.format(p))
#         pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#         col_mean = []
#         col_mle = []
#         for aa,protein in zip(*[pep.index,pep['Proteins']]):
#             if pd.notna(protein):
#                 u = protein.split(',')[0].lstrip('>')
#                 mean = uid_2_ts_mean[u]
#                 mle = uid_2_ts_mle[u]
#             else:
#                 mean,mle = pd.NA, pd.NA
#             col_mean.append(mean)
#             col_mle.append(mle)
#         pep['tumor_specificity_mean'] = col_mean
#         pep['tumor_specificity_mle'] = col_mle
#         pep.to_excel(writer,sheet_name=p)
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
#     # 18 neoantigens table
#     df = pd.read_csv('MS_common_occurence.txt',sep='\t',index_col=0)
#     df['uid'] = [aa_2_uid[item] for item in df.index]
#     df['tumor_specificity_mean'] = [uid_2_ts_mean[item] for item in df['uid']]
#     df['tumor_specificity_mle'] = [uid_2_ts_mle[item] for item in df['uid']]
#     df.to_excel(writer,sheet_name='18_common_neoantigens')
#     # add description
#     dic = {
#         'Sequence':'The amino acid sequence of the identified peptide',
#         'N-term cleavage window':'Sequence window from -15 to 15 around the N-terminal cleavage site of this peptide',
#         'C-term cleavage window':'Sequence window from -15 to 15 around the C-terminal cleavage site of this peptide',
#         'Amino acid before':'The amino acid in the protein sequence before the peptide',
#         'First amino acid':'The amino acid in the first position of the peptide sequence',
#         'Second amino acid':'The amino acid in the second position of the peptide sequence',
#         'Second last amino acid':'The amino acid in the second last position of the peptide sequence',
#         'Last amino acid':'The amino acid in the last position of the peptide sequence',
#         'Amino acid after':'The amino acid in the protein sequence after the peptide',
#         'N Count':'The number of instances of the "N" amino acid contained within the sequence, N indicates amino acid letter',
#         'Length':'The length of the sequence stored in the column "Sequence"',
#         'Missed cleavages':'Number of missed enzymatic cleavages',
#         'Mass':'Monoisotopic mass of the peptide',
#         'Proteins':'Identifiers of proteins this peptide is associated with',
#         'Leading razor protein':'Identifier of the leading protein in the protein group which uses this peptide for quantification. (Either unique or razor)',
#         'Start position':'Position of the first amino acid of this peptide in the protein sequence. (one-based)',
#         'End position':'Position of the last amino acid of this peptide in the protein sequence. (one-based)',
#         'Unique (Groups)':'When marked with "+", this particular peptide is unique to a single protein group in the proteinGroups file',
#         'Unique (Proteins)':'When marked with "+", this particular peptide is unique to a single protein sequence in the fasta file(s)',
#         'Charges':'All charge states that have been observed',
#         'PEP':'Posterior Error Probability of the identification. This value essentially operates as a p-value, where smaller is more significant',
#         'Score':'Highest Andromeda score for the associated MS/MS spectra',
#         'Experient [n]': 'Number of evidence entries for this "Experiment [n]"',
#         'Intensity':'Summed up eXtracted Ion Current (XIC) of all isotopic clusters associated with the identified AA sequence. In case of a labeled experiment this is the total intensity of all the isotopic patterns in the label cluster',
#         'Reverse':'When marked with "+", this particular peptide was found to be part of a protein derived from the reversed part of the decoy database. These should be removed for further data analysis',
#         'Potential contaminant':'When marked with '+', this particular peptide was found to be part of a commonly occurring contaminant. These should be removed for further data analysis',
#         'id':'A unique (consecutive) identifier for each row in the peptides table, which is used to cross-link the information in this table with the information stored in the other tables',
#         'Protein group IDs':'The identifiers of the protein groups this peptide was linked to, referenced against the proteinGroups table',
#         'Mod. peptide IDs':'Identifier(s) for peptide sequence(s), associated with the peptide, referenced against the corresponding modified peptides table',
#         'Evidence IDs':'Identifier(s) for analyzed peptide evidence associated with the protein group referenced against the evidence table',
#         'MS/MS IDs':'The identifiers of the MS/MS scans identifying this peptide, referenced against the msms table',
#         'Best MS/MS':'The identifier of the best (in terms of quality) MS/MS scan identifying this peptide, referenced against the msms table',
#         'Oxidation (M) site IDs':'Identifier(s) for site(s) associated with the protein group, which show(s) evidence of the modification, referenced against the appropriate modification site file',
#         'Taxonomy IDs':'Taxonomy identifiers',
#         'MS/MS Count':'The number of MS/MS evidence',
#         'tumor_specificity_mean':'The average raw read counts for each parental NeoJunction',
#         'tumor_specificity_mle':'The Tumor specificity score derived using Maximum Likelihood Estimation(MLE)',
#     }
#     pd.Series(data=dic,name='description').to_excel(writer,sheet_name='description')








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



'''neoantigen analysis'''
# 1. physicalchemical properties and occurence frequency
# snaf.downstream.analyze_neoantigens(freq_path='result/frequency_stage3_verbosity1_uid.txt',junction_path='result/burden_stage0.txt',total_samples=472,outdir='result',mers=None,fasta=False)
# snaf.run_dash_T_antigen(input_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/result/shared_vs_unique_neoantigen_all.txt',
#                         output_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/result')
# snaf.downstream.plot_umap_neoantigen(df_path='result/mer9_umap_embed_df.txt',outdir='result')
        




























