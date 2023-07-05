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


# # # get reduced junction
# df = snaf.get_reduced_junction_matrix(pc='counts.TCGA-SKCM.txt',pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')
# # df.to_csv('../ts/logic_gate/tcga_melanoma_df.txt',sep='\t')

# # run SNAF
# netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
# db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
# tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
# gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
# add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

# snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
# surface.initialize(db_dir=db_dir)

# show PMEL ENSG00000185664:E11.9-E12.2_55956189
# uid = 'ENSG00000185664:E11.9-E12.2_55956189'
# snaf.gtex_visual_combine(uid=uid,norm=False,outdir='Frank_inspection',tumor=df,group_by_tissue=False)
# snaf.gtex_visual_combine(uid=uid,norm=True,outdir='Frank_inspection',tumor=df,group_by_tissue=False)
# adata_gene = ad.read_h5ad('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ts/logic_gate/coding.h5ad')
# snaf.gtex_viewer.gtex_viewer_configuration(adata_gene)
# uid = 'ENSG00000185664'
# snaf.gtex_visual_combine(uid=uid,norm=False,outdir='Frank_inspection',tumor=df,group_by_tissue=False)
# snaf.gtex_visual_combine(uid=uid,norm=True,outdir='Frank_inspection',tumor=df,group_by_tissue=False)
# sys.exit('stop')

# 4 common neoantigens for immuno assay
# cand = pd.read_csv('result_new/T_candidates/T_antigen_candidates_all.txt',sep='\t',index_col=0)
# outdir = '4_common_inspect'
# uid = 'ENSG00000164175:E3.2_33963931-E4.2'
# criterion=[('netMHCpan_el', 0, '<=', 2),('deepimmuno_immunogenicity',1,'==','True')]
# criterion = [('netMHCpan_el', 0, '<=', 2)]
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

# # some stats
# burden = pd.read_csv('result_new/burden_stage3.txt',sep='\t',index_col=0)
# a = burden.iloc[-1,:-1]  
# print(a.mean(),a.max(),a.min())
'''
stage 0: 527.9216101694915 1549.0 28.0
stage 2: 1090.4300847457628 2981.0 75.0
stage 3: 915.2754237288135 2486.0 74.0
'''

# build a table with each sample, and the associated three stage burden, and the classification
lis = []
for s in [0,2,3]:
    burden = pd.read_csv('result_new/burden_stage{}.txt'.format(s),sep='\t',index_col=0)
    a = burden.iloc[-1,:-1]
    a.name = 'stage{}_burden'.format(s)
    lis.append(a)
stat_df = pd.concat(lis,axis=1)
classify = pd.read_csv('result_new/survival/groups.txt',sep='\t',index_col=0,header=None)
dic = classify[2].to_dict()
stat_df['identity'] = stat_df.index.map(dic).values
stat_df.to_csv('stat_table_number_of_burden.txt',sep='\t')
sys.exit('stop')
    

# do survival and mutation analysis
# survival = pd.read_csv('TCGA-SKCM.survival.tsv',sep='\t',index_col=0)
# for stage in [0,2,3]:
#     burden_df = pd.read_csv('result_new/burden_stage{}.txt'.format(stage),sep='\t',index_col=0)
#     burden_df_subset = burden_df.loc[burden_df.iloc[:,-1]!=0,:]
#     unique_n_neojunctions = list(set(burden_df_subset.index.tolist())) # 16800 9541 9323
#     burden_df.rename(columns=lambda x:'-'.join(x.split('-')[:4]),inplace=True)
#     burden_output,_ = snaf.survival_analysis(burden_df,survival,2,stratification_plot='result_new/survival/stage{}_stratify.pdf'.format(stage),
#                                              survival_plot='result_new/survival/stage{}_survival.pdf'.format(stage))
#     burden_output.to_csv('result_new/survival/burden{}_patient_high_low_group.txt'.format(stage),sep='\t')


# mutation = pd.read_csv('TCGA-SKCM.mutect2_snv.tsv',sep='\t',index_col=0)
# # burden3 = pd.read_csv('result_new/burden_stage3.txt',sep='\t',index_col=0)
# # burden3.rename(columns=lambda x:'-'.join(x.split('-')[:4]),inplace=True)
# # snaf.mutation_analysis(mode='compute',burden=burden3,mutation=mutation,output='result_new/survival/mutation.txt')
# # snaf.mutation_analysis(mode='plot',burden=burden3,mutation=mutation,output='result_new/survival/CAMKK2_mutation.txt',genes_to_plot=['CAMKK2'])
# mutation_camkk2 = mutation.loc[mutation['gene']=='CAMKK2',:].index.tolist()
# burden3 = pd.read_csv('result_new/survival/burden3_patient_high_low_group.txt',sep='\t',index_col=0)
# burden3_high = burden3.loc[burden3['identity']=='high',:]
# high_burden_samples = burden3.index.tolist()
# occur_in_high = list(set(mutation_camkk2).intersection(set(high_burden_samples)))
# print(len(occur_in_high))  # 18


# snaf.downstream.survival_regression(freq='result_new/frequency_stage3_verbosity1_uid_gene_symbol_coord_mean_mle.txt',remove_quote=True,
#                                     rename_func=lambda x:'-'.join(x.split('-')[:4]),survival='TCGA-SKCM.survival.tsv',
#                                     pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',outdir='result_new/survival',mode='binary')

# do DEG and GO
# snaf.downstream.prepare_DEG_analysis('result_new/burden_stage3.txt','result_new/survival/burden3_patient_high_low_group.txt',
#                                      rename_func=lambda x:'-'.join(x.split('-')[:4]),
#                                      outdir='result_new/survival',encoding={'low':'2','high':'1'})
# snaf.downstream.visualize_DEG_result('result_new/survival/DEGs-LogFold_0.0_adjp/GE.low_vs_high_mod.txt',up_cutoff=0.58,down_cutoff=-0.58,
#                                      mode='static',outdir='result_new/survival',genes_to_highlight=['LST1','HCST','IL32','CD3D','S100A8','MZB1','IGLC4','ADAM10','ARFGEF2','MIB1','KIF3B','TNPO1','PTPN11','ANKRD52','TGFBR1'])

# visualize RBP specifically
# rbps = ['EIF4G2', 'HLTF', 'DDX3X', 'NOLC1', 'G3BP2', 'XRN2', 'FAM120A', 'WDR3', 'NAA15', 'PNPT1', 'DDX21', 'XPO1', 'ZC3H11A', 'PUM1', 'ADAR', 'WDR43', 'EIF3A', 'PUM2', 'UCHL5', 'ZRANB2']
# snaf.downstream.visualize_DEG_result('result_new/survival/DEGs-LogFold_0.0_adjp/GE.low_vs_high_mod.txt',up_cutoff=0.58,down_cutoff=-0.58,
#                                      mode='static',outdir='result_new/survival',genes_to_highlight=rbps)

# snaf.downstream.prepare_GO_analysis('result_new/survival/DEGs-LogFold_0.0_adjp/GE.low_vs_high.txt',outdir='result_new/survival',direction='>',lc_cutoff=0.58,adjp_cutoff=0.05)
# snaf.downstream.visualize_GO_result(path_list=['result_new/survival/GO_Elite_result_BioMarkers/GO-Elite_results/CompleteResults/ORA/archived-20230528-192111/gene_list_up_in_low-BioMarkers.txt','result_new/survival/GO_Elite_result_GeneOntology/GO-Elite_results/CompleteResults/ORA/archived-20230528-192137/gene_list_up_in_low-GO.txt'],
#                                     skiprows_list=[16,17],category_list=['Gene-Set Name','Ontology Name'],outdir='result_new/survival',
#                                     mode='static',ontology_to_highlight={'Adult Peripheral Blood Activated T cell (PMID32214235 top 100)':'T cells','antigen binding':'antigen binding','complement activation':'Complement Activation','immune response':'immune response','humoral immune response':'humoral immune response'},ylims=(10e-85,10e-1))
# snaf.downstream.visualize_GO_result(path_list=['result_new/survival/GO_Elite_result_GeneOntology/GO-Elite_results/CompleteResults/ORA/archived-20230528-192717/gene_list_up_in_high-GO.txt'],
#                                     skiprows_list=[17],category_list=['Ontology Name'],outdir='result_new/survival',
#                                     mode='static',ontology_to_highlight={'cellular protein metabolic processs':'metabolic process','regulation of DNA endoreduplication':'DNA endoreduplication','cadherin binding':'cadherin binding','synapse part':'synapse part','cellular component organization or biogenesis':'cellular biogenesis','nucleoside-triphosphatase activity':'nucleoside-triphosphatase activity','positive regulation of macromolecule metabolic process':'regulation of metabolic'},ylims=(10e-35,10e-0))

# # coverage analysis
# snaf.downstream.get_coverage(t_result='result_new/T_candidates/T_antigen_candidates_all.txt',allele='A')

# common neoantigen characteristis
# freq3 = pd.read_csv('result_new/frequency_stage3_verbosity1_uid_gene_symbol_coord_mean_mle.txt',sep='\t',index_col=0)
# freq3_shared = freq3.loc[freq3['n_sample']>472*0.15,:]
# genes = set(freq3_shared['symbol'].values)
# with open('result_new/common/gene_list.txt','w') as f:
#     for gene in genes:
#         f.write('{}\n'.format(gene))
# snaf.downstream.visualize_GO_result(path_list=['result_new/common/GO_Elite_result_GeneOntology/GO-Elite_results/CompleteResults/ORA/archived-20230428-094602/gene_list-GO.txt'],
#                                     skiprows_list=[17],category_list=['Ontology Name'],outdir='result_new/common',
#                                     mode='static',ontology_to_highlight={'melanocyte differentiation':'melanocyte differentiation',
#                                     'melanosome membrane':'melanosome membrane','pigment cell differentiation':'pigment cell differentiation','melanin metabolic process':'melanin metabolic process'},ylims=(10e-5,10e-1),xlims=(-1,16))


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

# some statistics of identified hits
sr_str3_false = pd.read_csv('result_new/surface/B_candidates/sr_str3_report_None_False.txt',sep='\t',index_col=0)
false_uid = [j+','+e for j,e in zip(sr_str3_false['NeoJunction'],sr_str3_false['evidence'])]
all_uid_str3 = list(set(false_uid))   # 378

sr_str4_false = pd.read_csv('result_new/surface/B_candidates/sr_str4_report_None_False.txt',sep='\t',index_col=0)
false_uid = [j+','+e for j,e in zip(sr_str4_false['NeoJunction'],sr_str4_false['evidence'])]
all_uid_str4 = list(set(false_uid))   

str3_pass_str4 = list(set(all_uid_str3).intersection(set(all_uid_str4)))  # 37

sr_str5_false = pd.read_csv('result_new/surface/B_candidates/sr_str5_report_None_False.txt',sep='\t',index_col=0)
false_uid = [j+','+e for j,e in zip(sr_str5_false['NeoJunction'],sr_str5_false['evidence'])]
all_uid_str5 = list(set(false_uid))   

str4_pass_str5 = list(set(all_uid_str4).intersection(set(all_uid_str5)))  # 17

# show SIRPA example
# snaf.gtex_visual_combine(uid='ENSG00000198053:E7.2-E13.1_1915159',norm=False,outdir='result_new/surface',tumor=df,group_by_tissue=False)

# show EDA2R example
'''deletion and insertion can be wrong depending on which ref_seq you chose but by large the classification is correct, 
  extracellular prediction can be wrong if the topology is complicated, None in validation column in str4 may due to gtf issue
  the most inclusive result is None_False, if you need extracellular, go to None_True, if you need deletion or insertion, go to relevant subsets'''
# snaf.gtex_visual_combine(uid='ENSG00000131080:E6.1-E8.1',norm=False,outdir='result_new/surface',tumor=df,group_by_tissue=False)
# surface.run_dash_B_antigen(pkl='result_new/surface/surface_antigen_sr.p',candidates='result_new/surface/candidates_5_sr_None_False.txt',prediction_mode='short_read',
#                            python_executable='/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7')

lr_str3_false = pd.read_csv('result_new/surface/B_candidates/lr_str3_report_None_False.txt',sep='\t',index_col=0)
false_uid = [j+','+e for j,e in zip(lr_str3_false['NeoJunction'],lr_str3_false['evidence'])]
all_uid_str3 = list(set(false_uid))     # 1207

unique_valid_sr = set(sr_str5_false['peptide_sequence'].tolist())
unique_valid_lr = set(lr_str3_false['peptide_sequence'].tolist())
total_valid = unique_valid_sr.union(unique_valid_lr)  # 562

# show SLC45A2
# snaf.gtex_visual_combine(uid='ENSG00000164175:E3.2_33963931-E4.2',norm=False,outdir='result_new/surface',tumor=df,group_by_tissue=False)

# explore all other examples
'''when sorting by mean count in ascending order, common disquanlification reasons are
(1) although lowly in normal, not high in tumor either
(2) too disruptive sequence when considering intron retention or novel exon or other
(3) ref seq used internally were not comprehensive enough, so actually a documented isoform
(4) weired one like circular rna'''
# surface.run_dash_B_antigen(pkl='result_new/surface/surface_antigen_lr.p',candidates='result_new/surface/candidates_3_lr_None_False.txt',prediction_mode='long_read',
#                            python_executable='/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7')

# show four examples, ANO10, IGSF11, NALCN, MET_new
snaf.gtex_visual_combine(uid='ENSG00000092421:E22.1-E24.1_116468915',norm=False,outdir='result_new/surface',tumor=df,group_by_tissue=False)
snaf.gtex_visual_combine(uid='ENSG00000057019:E5.1-I5.1_98867438',norm=False,outdir='result_new/surface',tumor=df,group_by_tissue=False)
sys.exit('stop')
snaf.gtex_visual_combine(uid='ENSG00000160746:I21.1_43457363-E22.1',norm=False,outdir='result_new/surface',tumor=df,group_by_tissue=False)
snaf.gtex_visual_combine(uid='ENSG00000144847:E4.4-E8.1',norm=False,outdir='result_new/surface',tumor=df,group_by_tissue=False)
snaf.gtex_visual_combine(uid='ENSG00000102452:U0.1_101417117-E2.1',norm=False,outdir='result_new/surface',tumor=df,group_by_tissue=False)
snaf.gtex_visual_combine(uid='ENSG00000105976:E15.2-E17.2_116763110',norm=False,outdir='result_new/surface',tumor=df,group_by_tissue=False)
sys.exit('stop')


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
# ax.set_xlim([-2,12])
# sns.kdeplot(data=series_common.to_frame(),x='occurence',ax=ax,clip=[-2,15])
# sns.kdeplot(data=series_unique.to_frame(),x='occurence',ax=ax,clip=[-2,15])
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

# common versus unique event type
# freq = pd.read_csv('result_new/frequency_stage3_verbosity1_uid_gene_symbol_coord_mean_mle.txt',sep='\t',index_col=0)
# common_uid = list(set([item.split(',')[1] for item in freq.loc[freq['n_sample']>472*0.15,:].index]))
# unique_uid = list(set([item.split(',')[1] for item in freq.loc[freq['n_sample']==1,:].index]))
# snaf.downstream.plot_event_type(pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',uids={'common':common_uid,'unique':unique_uid},rel=True,outdir='result_new')

'''
inspect those candidates for junction validity and tumor specificity
'''
control_bam_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A208-11A-51R-A157-07.bam',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-Kidney/KIRC/TCGA-CJ-5680-11A-01R-1541-07.bam',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-STAD/TCGA-BR-6453-11A-01R-1802-13.bam',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-ESCA/TCGA-L5-A4OG-11A-12R-A260-31.bam',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-THCA/TCGA-EL-A3ZP-11A.bam',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A0H7-11A-13R-A089-07.bam']
control_bai_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A208-11A-51R-A157-07.bam.bai',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-Kidney/KIRC/TCGA-CJ-5680-11A-01R-1541-07.bam.bai',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-STAD/TCGA-BR-6453-11A-01R-1802-13.bam.bai',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-ESCA/TCGA-L5-A4OG-11A-12R-A260-31.bam.bai',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-THCA/TCGA-EL-A3ZP-11A.bam.bai',
                 '/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A0H7-11A-13R-A089-07.bam.bai']

# control_bam_path_list = ['/data/salomonis-archive/BAMs/PublicDatasets/E-MTAB-2836-Grch38_Deep-Healthy-PanTissue/ERR315409_1.bam']
# control_bai_path_list = ['/data/salomonis-archive/BAMs/PublicDatasets/E-MTAB-2836-Grch38_Deep-Healthy-PanTissue/ERR315409_1.bam.bai']

# control_bam_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A208-11A-51R-A157-07.bam']
# control_bai_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-BREAST-CANCER/TCGA-files-Ens91/bams/TCGA-BH-A208-11A-51R-A157-07.bam.bai']

def unit_run(uid,base_sample,region,task_name):
    uid = uid
    sample = base_sample + '.bed'
    region = region
    bam_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-SKCM/TCGA_SKCM-BAMs/bams/{}.bam'.format(base_sample)] + control_bam_path_list
    bai_path_list = ['/data/salomonis-archive/BAMs/NCI-R01/TCGA/TCGA-SKCM/TCGA_SKCM-BAMs/bams/{}.bam.bai'.format(base_sample)] + control_bai_path_list
    sif_anno_path = '/data/salomonis2/software/ggsashimi'
    outdir = 'result_new/common/sashimi'
    bam_contig_rename = [False] + [False] * 5
    criterion=[('netMHCpan_el', 0, '<=', 2),('deepimmuno_immunogenicity',1,'==','True')]

    snaf.gtex_visual_combine_plotly(uid=uid,outdir='result_new/common',norm=False,tumor=df)
    snaf.gtex_visual_combine_plotly(uid=uid,outdir='result_new/common',norm=True,tumor=df)
    snaf.JunctionCountMatrixQuery.deserialize('result_new/after_prediction.p').visualize(uid=uid,sample=sample,outdir='result_new/common',criterion=criterion)

    # dff = snaf.gtex_visual_combine(uid=uid,outdir='Frank_inspection',norm=False,tumor=df)
    # dff.to_csv('HAAASFETL_adata_df.txt',sep='\t')
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

'''find sashimi plots for those common neoantigens'''
freq3 = pd.read_csv('result_new/frequency_stage3_verbosity1_uid_gene_symbol_coord_mean_mle.txt',sep='\t',index_col=0)
freq3 = freq3.loc[freq3['n_sample']>472*0.15,:]
occurence = pd.read_csv('MS_common_occurence.txt',sep='\t',index_col=0)
occurence = occurence.loc[occurence['occurence']>0,:].index.tolist()
freq3['peptide'] = [item.split(',')[0] for item in freq3.index]
freq3['junction'] = [item.split(',')[1] for item in freq3.index]
freq3['ms_evidence'] = freq3['peptide'].isin(set(occurence))
freq3 = freq3.loc[freq3['ms_evidence'],:]
col = []
for item in freq3['symbol']:
    if '-AS' in item or 'LINC' in item or 'unknown_gene' in item or 'LOC' in item:
        col.append(False)
    else:
        col.append(True)
freq3 = freq3.loc[col,:]
freq3 = freq3.loc[freq3['tumor_specificity_mean']<1,:]
freq3.drop_duplicates(subset=['junction'],inplace=True)

count = snaf.remove_trailing_coord('counts.TCGA-SKCM.txt',sep='\t')
col = []
from ast import literal_eval
freq3['samples'] = [literal_eval(item) for item in freq3['samples']]
for row in freq3.itertuples():
    uid = row.junction
    base_sample = row.samples[0].split('.')[0]
    c = count.at[uid,base_sample+'.bed']
    col.append(c)
freq3['count'] = col
# freq3.to_csv('result_new/common/sashimi_checklist.txt',sep='\t')

for row in freq3.itertuples():
    uid = row.junction
    base_sample = row.samples[0].split('.')[0]
    unit_run(uid,base_sample,None,None)
sys.exit('stop')






# specifically, generate publication-quality sashimi for 8 validated hits
'''
VAPGEAKNL
ENSG00000185989:E5.1_114128541-E18.1
Mel4
TCGA-DA-A1HY-06A-11R-A18T-07
RASA3
gene is highly expressed, but splicing is unique, very tumor specific
adipose tissue ERR315343 use this as control (not using) instead using breast tcga paratumor TCGA-BH-A0H7-11A-13R-A089-07
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
# unit_run(uid='ENSG00000132185:E2.3-E2.5',
#          base_sample='TCGA-EB-A3HV-01A-11R-A21D-07',
#          region=flank_chrom('chr1:161710492-161710760',(5000,1000)),
#          task_name='_HAAASFETL')

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
# unit_run(uid='ENSG00000100225:E8.2-E9.2_32491143',
#          base_sample='TCGA-ER-A2NE-06A-21R-A18T-07',
#          region=flank_chrom('chr22:32487828-32491143',(500,500)),
#          task_name='_KEKLDQLVY')


'''
additionally, few needs to be visualized as they can predict survival
ENSG00000137265:E4.2-E6.3
TCGA-D9-A149-06A-11R-A18S-07
IRF4
chr6:395935-397111
seems to be gene-driven
'''
# unit_run(uid='ENSG00000137265:E4.2-E6.3',
#          base_sample='TCGA-D9-A149-06A-11R-A18S-07',
#          region=flank_chrom('chr6:395935-397111',(1000,1000)),
#          task_name='_IRF4')

'''
additionally, few needs to be visualized as they can predict survival
ENSG00000164175:E3.2_33963931-E4.2
TCGA-D9-A148-06A-11R-A18S-07
SLC45A2
chr5:33954504-33963931
'''
# unit_run(uid='ENSG00000164175:E3.2_33963931-E4.2',
#          base_sample='TCGA-D9-A148-06A-11R-A18S-07',
#          region=flank_chrom('chr5:33954504-33963931',(1000,1000)),
#          task_name='_SLC45A2')



'''compare against the proliferative tissues'''
# # segragate proliferating and non-proliferating samples 
# count = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/proliferating/counts.original.txt',sep='\t',index_col=0)
# inventory = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/proliferating/inventory.txt',sep='\t',index_col=0)
# proliferating = [item + '_secondAligned.sortedByCoord.out.bed' for item in inventory.loc[inventory['cond'],:].index.tolist()]
# non_proliferating = [item + '_secondAligned.sortedByCoord.out.bed' for item in inventory.loc[~inventory['cond'],:].index.tolist()]
# proliferating.remove('SRR18143890_secondAligned.sortedByCoord.out.bed')  # not successfully aligned
# count_pro = count.loc[:,proliferating]   # [822922 rows x 125 columns]
# count_pro = count_pro.loc[count_pro.sum(axis=1)!=0,:]  # [690742 rows x 125 columns]
# count_non = count.loc[:,non_proliferating]  # [822922 rows x 44 columns]
# count_non = count_non.loc[count_non.sum(axis=1)!=0,:]  # [676570 rows x 44 columns]
# count_pro.to_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/proliferating/count_pro.txt',sep='\t')
# count_non.to_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/proliferating/count_non.txt',sep='\t')

# get reduced junction
df = snaf.get_reduced_junction_matrix(pc='counts.TCGA-SKCM.txt',pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')

# # run SNAF but just get the maxmin text file
# # /data/salomonis-archive/FASTQs/PublicDatasets/Bulk-RNASeq/ProliferatingSkinCell
# netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
# db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
# tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
# gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
# proliferate_db = snaf.remove_trailing_coord('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/proliferating/count_pro.txt')
# proliferate_db.rename(columns=lambda x:x+'_proliferate',inplace=True)
# non_proliferate_db = snaf.remove_trailing_coord('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/proliferating/count_non.txt')
# non_proliferate_db.rename(columns=lambda x:x+'_non_proliferate',inplace=True)
# add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db,'proliferate':proliferate_db,'non_proliferate':non_proliferate_db}

# snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
# surface.initialize(db_dir=db_dir)
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=30,add_control=add_control,outdir='result_proliferate_new',filter_mode='maxmin')


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
reduce valid Neojunction from 16799 to 15017 because they are present in added control proliferate
reduce valid Neojunction from 15017 to 14553 because they are present in added control non_proliferate
'''

# get gene_list that are removed in last step
# filter_df = pd.read_csv('result_proliferate_new/NeoJunction_statistics_maxmin.txt',sep='\t',index_col=0)
# filter_df_removed = filter_df.loc[(filter_df['cond']==True)&(filter_df['cond_add_tcga_control']==True)&(filter_df['cond_add_gtex_skin']==True)&(filter_df['cond_add_proliferate']==False),:]
# filter_df_removed.to_csv('result_proliferate_new/filter_df_filtered.txt',sep='\t')
# sns.histplot(filter_df_removed['mean_add_proliferate'].values)
# plt.savefig('result_proliferate_new/new_mean_proliferate.pdf',bbox_inches='tight')
# plt.close()
# genes = list(set([item.split(':')[0] for item in filter_df_removed.index]))
# symbols = snaf.downstream.ensemblgene_to_symbol(query=genes,species='human')
# with open('result_proliferate_new/gene_list.txt','w') as f:
#     for s in symbols:
#         if s != 'unknown_gene':
#             f.write('{}\n'.format(s))
# snaf.visualize_GO_result(path_list=['/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_proliferate_new/GO_Elite_result_GeneOntology/GO-Elite_results/CompleteResults/ORA/archived-20230423-141705/gene_list-GO.txt'],
#                          skiprows_list=[17],category_list=['Ontology Name'],outdir='result_proliferate_new',mode='static',ylims=(10e-25,10e-1),
#                          ontology_to_highlight=
#                              {'cell cycle':'cell cycle',
#                               'cell division':'cell division',
#                               'cell cycle process':'cell cycle process',
#                               'chromosome segregation':'chromosome segregation',
#                               'spindle':'spindle',
#                               'cytoskeleton':'cytoskeleton',
#                               'regulation of cell cycle process':'regulation of cell cycle process'})





'''write peptide.txt to a merged xlsx file for publication, supp3 table'''
# os.chdir('/data/salomonis-archive/MS/melanoma/raw')
# all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
# all_patients.pop(all_patients.index('Mel-16'))
# df = pd.read_csv('result_new/frequency_stage3_verbosity1_uid_gene_symbol_coord_mean_mle.txt',sep='\t',index_col=0)
# df['uid'] = [item.split(',')[1] for item in df.index]
# df['aa'] = [item.split(',')[0] for item in df.index]
# uid_2_ts_mean = pd.Series(data=df['tumor_specificity_mean'].values,index=df['uid'].values).to_dict()
# uid_2_ts_mle = pd.Series(data=df['tumor_specificity_mle'].values,index=df['uid'].values).to_dict()
# aa_2_uid = pd.Series(data=df['uid'].values,index=df['aa'].values).to_dict()
# with pd.ExcelWriter('supp_table_ms.xlsx') as writer:
#     for p in all_patients:
#         os.chdir('/data/salomonis-archive/MS/melanoma/raw/{}/combined/txt'.format(p))
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
#     os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma')
#     # common neoantigens table
#     df = pd.read_csv('MS_common_occurence.txt',sep='\t',index_col=0)
#     df['uid'] = [aa_2_uid[item] for item in df.index]
#     df['tumor_specificity_mean'] = [uid_2_ts_mean[item] for item in df['uid']]
#     df['tumor_specificity_mle'] = [uid_2_ts_mle[item] for item in df['uid']]
#     df.to_excel(writer,sheet_name='common_neoantigens')
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



'''neoantigen analysis'''
# snaf.downstream.analyze_neoantigens(freq_path='result_new/frequency_stage3_verbosity1_uid.txt',junction_path='result_new/burden_stage0.txt',total_samples=472,outdir='result_new',mers=None,fasta=False)
# snaf.run_dash_T_antigen(input_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new/shared_vs_unique_neoantigen_all.txt',
#                         output_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/TCGA_melanoma/result_new')
snaf.downstream.plot_umap_neoantigen(df_path='result_new/mer9_umap_embed_df.txt',outdir='result_new')
        




























