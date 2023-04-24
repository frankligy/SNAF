#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
sys.path.insert(0,'/data/salomonis2/software')
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import re
import anndata as ad

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# the following code is trying to obtain a dictionary between srr_id and the patient_id and their HLA allele present in supplemental from the paper
rna = pd.read_csv('srr_list.txt',sep='\t',index_col=0,header=None)
pat = re.compile(r'Pat\d{2,3}')
patient = []
for item in rna.index:
    if not re.search(pat,item):
        raise Exception('check,stop')
    else:
        patient.append(re.search(pat,item).group(0))
rna['patient'] = patient
rna.set_index(keys='patient',inplace=True)
rna['full_name'] = [item + '_secondAligned.sortedByCoord.out.bed' for item in rna[1]]  
patient_dict = rna['full_name'].to_dict()   # {Pat40:SRR2660032_secondAligned.sortedByCoord.out.bed}

hla = pd.read_csv('all_patients_hla.txt',sep='\t',index_col=0).loc[:,['hla.a1','hla.a2','hla.b1','hla.b2','hla.c1','hla.c2']]
hla = hla.loc[hla.index.isin(patient_dict.keys()),:]
hla_dict = {}      # {Pat40:hla_string}
for row in hla.iterrows():
    hla_dict[row[0]] = ','.join(snaf.snaf.hla_formatting(row[1].tolist(),pre_type='netMHCpan_input',post_type='netMHCpan_output'))

full_dict = {}   # {SRR2660032_secondAligned.sortedByCoord.out.bed:hla_string}
for p,f in patient_dict.items():
    full_dict[f] = hla_dict[p]

# # run SNAF
# df = snaf.get_reduced_junction_matrix(pc='counts.original.txt',pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')
# netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
# db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
# tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
# gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
# add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

# snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
# surface.initialize(db_dir=db_dir)

# # T antigen
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=30,add_control=add_control,outdir='result_new',filter_mode='maxmin')
# hlas = [hla_string.split(',') for hla_string in df.columns.map(full_dict)]
# jcmq.run(hlas=hlas,outdir='./result_new')
# snaf.JunctionCountMatrixQuery.generate_results(path='./result_new/after_prediction.p',outdir='./result_new')

# # survival analysis
# survival = pd.read_csv('all_patients_hla.txt',sep='\t',index_col=0) 
reverse_patient_dict = {v:k for k,v in patient_dict.items()}
# for stage in [0,2,3]: 
#     burden = pd.read_csv('result_new/burden_stage{}.txt'.format(stage),sep='\t',index_col=0)
#     burden_subset = burden.loc[burden.iloc[:,-1]!=0,:]
#     unique_n_neojunctions = list(set(burden_subset.index.tolist())) # 15000, 9282, 9089
#     burden.columns = burden.columns.map(reverse_patient_dict).values
#     burden_output,quantiles = snaf.survival_analysis(burden,survival,n=2,stratification_plot='result_new/survival/stage{}_stratify.pdf'.format(stage),survival_plot='result_new/survival/stage{}_survival.pdf'.format(stage),survival_duration='overall_survival',survival_event='dead')
#     burden_output.to_csv('result_new/survival/burden{}_patient_high_low_group.txt'.format(stage),sep='\t')


# snaf.downstream.survival_regression(freq='result_new/frequency_stage3_verbosity1_uid_gene_symbol_coord_mean_mle.txt',remove_quote=True,
#                                     rename_func=lambda x:reverse_patient_dict[x],survival='all_patients_hla.txt',survival_duration='overall_survival',survival_event='dead',
#                                     pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt',outdir='result_new/survival',mode='binary')

df_tcga = pd.read_csv('../TCGA_melanoma/result_new/survival/survival_regression_binary_final_results.txt',sep='\t',index_col=0)
df_van = pd.read_csv('result_new/survival/survival_regression_binary_final_results.txt',sep='\t',index_col=0)
cross_df = df_tcga.join(df_van,how='inner',lsuffix='_tcga',rsuffix='_van')
cross_df.to_csv('result_new/survival/cross_df.txt',sep='\t')

sys.exit('stop')


# 2. embed neoantigen to umap
# snaf.downstream.analyze_neoantigens(freq_path='result/frequency_stage3_verbosity1_uid.txt',junction_path='result/burden_stage0.txt',total_samples=39,outdir='result',mers=None,fasta=False)
# snaf.run_dash_T_antigen(input_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/immunotherapy/allen/result/shared_vs_unique_neoantigen_all.txt',
#                         output_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/immunotherapy/allen/result')
# snaf.downstream.plot_umap_neoantigen(df_path='result/mer9_umap_embed_df.txt',outdir='result')


# survival correlation
# freq = 0.20
# tcga = pd.read_csv('survival_correlation/frequency_stage2_verbosity1_uid_TCGA_SKCM.txt',sep='\t',index_col=0)
# tcga.index = [item.replace(',','.').replace(':','.').replace('-','.') for item in tcga.index]
# tcga = tcga.loc[tcga['n_sample']>freq*472,:]
# tcga_dict = (tcga['n_sample']/472).to_dict()
# allen = pd.read_csv('survival_correlation/frequency_stage2_verbosity1_uid_van_allen.txt',sep='\t',index_col=0)
# allen.index = [item.replace(',','.').replace(':','.').replace('-','.') for item in allen.index]
# allen = allen.loc[allen['n_sample']>freq*39,:]
# allen_dict = (allen['n_sample']/39).to_dict()

# surv_corr_tcga = pd.read_csv('survival_correlation/SurvResult_dPSI_stage2_TCGA_SKCM.txt',sep='\t',index_col=0)
# surv_corr_tcga = surv_corr_tcga.loc[surv_corr_tcga.index.isin(tcga.index),:]
# surv_corr_allen = pd.read_csv('survival_correlation/SurvResult_dPSI_stage2_VanAllen.txt',sep='\t',index_col=0)
# surv_corr_allen = surv_corr_allen.loc[surv_corr_allen.index.isin(allen.index),:]

# surv_corr_merge = surv_corr_tcga.join(other=surv_corr_allen,how='inner',lsuffix='_tcga',rsuffix='_allen')
# surv_corr_merge = surv_corr_merge.loc[(surv_corr_merge['LRT_tcga']<0.05)&(surv_corr_merge['LRT_allen']<0.05),:]
# surv_corr_merge['freq_tcga'] = surv_corr_merge.index.map(tcga_dict).values
# surv_corr_merge['freq_allen'] = surv_corr_merge.index.map(allen_dict).values
# surv_corr_merge.to_csv('survival_correlation/freq_{}_result.txt'.format(freq),sep='\t')

# df = pd.read_csv('survival_correlation/freq_0.05_result.txt',sep='\t',index_col=0)
# fig,ax = plt.subplots()
# sp = ax.scatter(df['Zscore_tcga'].values,df['Zscore_allen'].values,s=20,c=df['freq_tcga'].values)
# for direc in ['top','right']:
#     ax.spines[direc].set_visible(False)
# for direc in ['bottom','left']:
#     ax.spines[direc].set_position('zero')
# ax.set_xlim((-6.01,6.01))
# ax.set_ylim((-5.01,5.01))
# ax.set_xticks(np.arange(-6,7,1))
# ax.set_xticklabels(np.arange(-6,7,1))
# ax.set_yticks(np.arange(-5,6,1))
# ax.set_yticklabels(np.arange(-5,6,1))
# ax.grid(alpha=0.5)
# plt.colorbar(mappable=sp)
# plt.savefig('survival_correlation/freq_0.05.pdf',bbox_inches='tight')
# plt.close()


# modify the supp2_table
with pd.ExcelWriter('/data/salomonis2/LabFiles/Frank-Li/neoantigen/immunotherapy/allen/survival_correlation/supp_table2.xlsx') as writer:
    # surv result tcga
    df = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/result/frequency_stage2_verbosity1_uid_gene_symbol_coord_mean_mle.txt',sep='\t',index_col=0)
    item_2_n = df['n_sample'].to_dict()
    item_2_ts_mean = df['tumor_specificity_mean'].to_dict()
    item_2_ts_mle = df['tumor_specificity_mle'].to_dict()
    sheet = pd.read_csv('survival_correlation/SurvResult_dPSI_stage2_TCGA_SKCM.txt',sep='\t',index_col=0)
    identifier = []
    col_n = []
    col_ts_mean = []
    col_ts_mle = []
    for index in tqdm(sheet.index,total=sheet.shape[0]):
        n_parts = len(index.split('.'))
        part1 = index.split('.')[0]
        part2 = index.split('.')[1]
        part3 = '.'.join(index.split('.')[2:4])
        if n_parts == 6:
            part4 = '.'.join(index.split('.')[4:])
        elif n_parts == 7:
            part4 = '.'.join(index.split('.')[4:])
            ensg = part4.split('.')[0]
            exons = '.'.join(part4.split('.')[1:])
            part4 = ':'.join([ensg,exons])
        item = part1 + ',' + part2 + ':' + part3 + '-' + part4
        n = item_2_n[item]
        ts_mean = item_2_ts_mean[item]
        ts_mle = item_2_ts_mle[item]
        identifier.append(item)
        col_n.append(n)
        col_ts_mean.append(ts_mean)
        col_ts_mle.append(ts_mle)
    sheet['identifier'] = identifier
    sheet['n_sample'] = col_n
    sheet['tumor_specificity_mean'] = col_ts_mean
    sheet['tumor_specificity_mle'] = col_ts_mle
    sheet.to_excel(writer,sheet_name='survival_result_tcga')
    # surv result van allen
    df = pd.read_csv('/data/salomonis2/LabFiles/Frank-Li/neoantigen/immunotherapy/allen/result/frequency_stage2_verbosity1_uid_gene_symbol_coord_mean_mle.txt',sep='\t',index_col=0)
    item_2_n = df['n_sample'].to_dict()
    item_2_ts_mean = df['tumor_specificity_mean'].to_dict()
    item_2_ts_mle = df['tumor_specificity_mle'].to_dict()
    sheet = pd.read_csv('survival_correlation/SurvResult_dPSI_stage2_VanAllen.txt',sep='\t',index_col=0)
    identifier = []
    col_n = []
    col_ts_mean = []
    col_ts_mle = []
    for index in tqdm(sheet.index,total=sheet.shape[0]):
        n_parts = len(index.split('.'))
        part1 = index.split('.')[0]
        part2 = index.split('.')[1]
        part3 = '.'.join(index.split('.')[2:4])
        if n_parts == 6:
            part4 = '.'.join(index.split('.')[4:])
        elif n_parts == 7:
            part4 = '.'.join(index.split('.')[4:])
            ensg = part4.split('.')[0]
            exons = '.'.join(part4.split('.')[1:])
            part4 = ':'.join([ensg,exons])
        item = part1 + ',' + part2 + ':' + part3 + '-' + part4
        n = item_2_n[item]
        ts_mean = item_2_ts_mean[item]
        ts_mle = item_2_ts_mle[item]
        identifier.append(item)
        col_n.append(n)
        col_ts_mean.append(ts_mean)
        col_ts_mle.append(ts_mle)
    sheet['identifier'] = identifier
    sheet['n_sample'] = col_n
    sheet['tumor_specificity_mean'] = col_ts_mean
    sheet['tumor_specificity_mle'] = col_ts_mle
    sheet.to_excel(writer,sheet_name='survival_result_van_allen')
    # biomarker
    freq = pd.read_csv('survival_correlation/freq_0.05_result.txt',sep='\t',index_col=0)
    freq.to_excel(writer,sheet_name='biomarkers')





























