import os
import sys
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import subprocess

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


'''
First, we need to request the permission on dbGAP to download the 472 TCGA SKCM patients' BAM files. 
Then we can run the first step of AltAnalyze, which has been detailed in the tutorial, after that, you will have a count matrix,
we provide this count matrix in synapse so that we can start the following analysis
The correponding input files and example output files are in https://www.synapse.org/#!Synapse:syn32057190
'''

df = pd.read_csv('./counts_matrix_tcga_melanoma.txt',sep='\t',index_col=0)

# run SNAF
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
db_dir = '/data/salomonis2/LabFiles/Frank-Li/refactor/data'

snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='netMHCpan',software_path=netMHCpan_path)
surface.initialize(db_dir=db_dir)

'''T cell neoantigen'''
jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=50)
sample_to_hla = pd.read_csv('../HLA_genotype/sample_hla.txt',sep='\t',index_col=0)['hla'].to_dict()
hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
jcmq.run(hlas=hlas,outdir='./result')
snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')




'''B cell neoantigen'''
membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df)
surface.run(membrane_tuples,outdir='result',tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm')
surface.generate_results(pickle_path='./result/surface_antigen.p',outdir='result',strigency=4,gtf='./SQANTI-all/collapse_isoforms_classification.filtered_lite.gtf')   # './SQANTI-all/collapse_isoforms_classification.filtered_lite.gtf'
surface.run_dash_B_antigen(pkl='result/surface_antigen.p',candidates='result/candidates_5.txt',python_executable='/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7')

print(snaf.uid_to_coord('ENSG00000164175:E3.2_33963931-E4.2'))
snaf.gtex_visual_combine('ENSG00000164175:E3.2_33963931-E4.2',norm=False,outdir='result',tumor=df)   # ENSG00000198053:E7.2-E13.1_1915159
snaf.gtex_visual_subplots('ENSG00000198053:E7.2-E13.1_1915159',norm=True,outdir='result')

'''T antigen MS validation'''
df = pd.read_csv('result/frequency_stage3_verbosity1_uid_gene_symbol.txt',sep='\t',index_col=0)
df_common = df.loc[df['n_sample']>282,:]
df_unique = df.loc[df['n_sample']==1,:]
with open('result/MS_validation_common.fasta','w') as f:
    for row in df_common.itertuples():
        peptide, uid = row.Index.split(',')
        f.write('>{}\n{}\n'.format(uid,peptide))
with open('result/MS_validation_unique.fasta','w') as f:
    for row in df_unique.itertuples():
        peptide, uid = row.Index.split(',')
        f.write('>{}\n{}\n'.format(uid,peptide))
for f in ['MS_validation_common','MS_validation_unique']:
    snaf.proteomics.remove_redundant('./result/{}.fasta'.format(f),'./result/{}_unique.fasta'.format(f))
    snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot_9_10_mers_unique.fasta',
                                      fa2_path='./result/{}_unique.fasta'.format(f),outdir='./result',write_unique2=True,prefix='{}_'.format(f))

db_common = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/fasta/MS_validation_common_unique2.fasta']
db_unique = ['/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/fasta/MS_validation_unique_unique2.fasta']

'''
common
'''
os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common')
all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
for p in all_patients:
    os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common/{}'.format(p))
    all_samples = subprocess.run('for file in *.raw; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    inputs = []
    pwd = os.getcwd()
    for s in all_samples:
        inputs.append(os.path.join(pwd,s))
    snaf.proteomics.set_maxquant_configuration(dbs=db_common,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                               outdir=pwd,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)


'''
unique
'''
os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_unique')
all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
for p in all_patients:
    os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_unique/{}'.format(p))
    all_samples = subprocess.run('for file in *.raw; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
    inputs = []
    pwd = os.getcwd()
    for s in all_samples:
        inputs.append(os.path.join(pwd,s))
    snaf.proteomics.set_maxquant_configuration(dbs=db_unique,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,
                                               outdir=pwd,minPepLen=8,minPeptideLengthForUnspecificSearch=8,maxPeptideLengthForUnspecificSearch=25)


'''
plot
'''
os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common')
all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
fig,ax = plt.subplots()
n_common = []
n_unique = []
from scipy.stats import mannwhitneyu
all_patients.pop(all_patients.index('Mel-16'))
for p in all_patients:
    os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common/{}/combined/txt'.format(p))
    pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
    pep = pep.loc[pep['Proteins'].notna(),:]
    n = int(pep.shape[0])/114
    n_common.append(n)
    os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_unique/{}/combined/txt'.format(p))
    pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
    pep = pep.loc[pep['Proteins'].notna(),:]
    n = int(pep.shape[0])/25188
    n_unique.append(n)  
stats = mannwhitneyu(n_common,n_unique)
print(stats)
number = len(all_patients)
ax.bar(x=[i for i in range(1,1+3*(number-1)+1,3)],height=n_common,label='common neoantigen')
ax.bar(x=[i for i in range(2,2+3*(number-1)+1,3)],height=n_unique,label='unique neoantigen')
xticks = []
for ic,iu in zip([i for i in range(1,1+3*(number-1)+1,3)],[i for i in range(2,2+3*(number-1)+1,3)]):
    xticks.append((ic+iu)/2)
ax.set_xticks(xticks)
ax.set_xticklabels(all_patients,fontsize=4,rotation=90)
ax.set_xlabel('patients',fontsize=6)
ax.set_ylabel('MS recovery rate',fontsize=6)
ax.legend(loc='upper left',bbox_to_anchor=(1,1),frameon=False)
os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
plt.savefig('MS_plot.pdf',bbox_inches='tight')
plt.close()

'''
occurence
'''
dict_common = {}
with open('result/MS_validation_common_unique2.fasta') as f:
    for line in f:
        if not line.startswith('>'):
            dict_common[line.rstrip('\n')] = 0
dict_unique = {}
with open('result/MS_validation_unique_unique2.fasta') as f:
    for line in f:
        if not line.startswith('>'):
            dict_unique[line.rstrip('\n')] = 0
os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common')
all_patients = subprocess.run('for file in *; do echo $file; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
all_patients.pop(all_patients.index('Mel-16'))
for p in all_patients:
    os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_common/{}/combined/txt'.format(p))
    pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
    pep = pep.loc[pep['Proteins'].notna(),:]
    for item in pep.index:
        try:
            dict_common[item] += 1
        except KeyError:
            continue
    os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/MS/raw_files_unique/{}/combined/txt'.format(p))
    pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
    pep = pep.loc[pep['Proteins'].notna(),:]
    for item in pep.index:
        try:
            dict_unique[item] += 1
        except KeyError:
            continue    
os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis')
series_common = pd.Series(dict_common)
series_common.name = 'occurence'
series_common.to_csv('MS_common_occurence.txt',sep='\t')
series_unique = pd.Series(dict_unique)
series_unique.name = 'occurence'
series_unique.to_csv('MS_unique_occurence.txt',sep='\t')
fig,ax = plt.subplots()
sns.ecdfplot(data=series_common.to_frame(),x='occurence',ax=ax)
sns.ecdfplot(data=series_unique.to_frame(),x='occurence',ax=ax)
import matplotlib.patches as mpatches
ax.legend(handles=[mpatches.Patch(color=i) for i in ['#4B71B0','#DE8353']],labels=['common neoantigen','unique neoantigen'],loc='upper left',bbox_to_anchor=(1,1),frameon=False)
#ax.text(x=0.5,y=0.1,s='Mann Whitney p-value: {}'.format(stats),transform=ax.transAxes)
plt.savefig('MS_occurence.pdf',bbox_inches='tight');plt.close()

'''common neoantigen GO analysis'''
df = pd.read_csv('result/toppfun_stage3_0.6.txt',sep='\t',encoding = 'unicode_escape')
df = df.loc[df['Category']=='GO: Biological Process',:]
df['convert'] = -np.log10(df['q-value FDR B&H'].values)
df = df.iloc[:5,:]
df.sort_values(by='convert',axis=0,ascending=True,inplace=True)
fig,ax = plt.subplots()
ax.barh(y=np.arange(df.shape[0]),width=[item for item in df['convert']],color='#FF9A91')
ax.set_yticks(np.arange(df.shape[0]))
ax.set_yticklabels([item for item in df['Name']])
ax.set_title('GO: Biological Process')
ax.set_xlabel('-Log10(adjusted_pval)')
plt.savefig('go_biological_process.pdf',bbox_inches='tight')
plt.close()

df = pd.read_csv('result/toppfun_stage3_0.6.txt',sep='\t',encoding = 'unicode_escape')
df = df.loc[df['Category']=='GO: Cellular Component',:]
df['convert'] = -np.log10(df['q-value FDR B&H'].values)
df = df.iloc[:5,:]
df.sort_values(by='convert',axis=0,ascending=True,inplace=True)
fig,ax = plt.subplots()
ax.barh(y=np.arange(df.shape[0]),width=[item for item in df['convert']],color='#FF9A91')
ax.set_yticks(np.arange(df.shape[0]))
ax.set_yticklabels([item for item in df['Name']])
ax.set_title('GO: Cellular Component')
ax.set_xlabel('-Log10(adjusted_pval)')
plt.savefig('go_cellular_component.pdf',bbox_inches='tight')
plt.close()

df = pd.read_csv('result/toppfun_stage3_0.6.txt',sep='\t',encoding = 'unicode_escape')
df = df.loc[df['Category']=='Pathway',:]
df['convert'] = -np.log10(df['q-value FDR B&H'].values)
df = df.iloc[:5,:]
df.sort_values(by='convert',axis=0,ascending=True,inplace=True)
fig,ax = plt.subplots()
ax.barh(y=np.arange(df.shape[0]),width=[item for item in df['convert']],color='#FF9A91')
ax.set_yticks(np.arange(df.shape[0]))
ax.set_yticklabels([item for item in df['Name']])
ax.set_title('Pathway')
ax.set_xlabel('-Log10(adjusted_pval)')
plt.savefig('pathway.pdf',bbox_inches='tight')
plt.close()


'''patient analysis'''
# 1. survival analysis
survival = pd.read_csv('TCGA-SKCM.survival.tsv',sep='\t',index_col=0)  # 463
burden = pd.read_csv('result/burden_stage0.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
burden_output,quantiles = snaf.survival_analysis(burden,survival,n=2,stratification_plot='result/stage0_stratify.pdf',survival_plot='result/stage0_survival.pdf')
burden_output.to_csv('result/to_nathan_stage0_neojunction_encode.txt',sep='\t')


# 2. mutation analysis
mutation = pd.read_csv('TCGA-SKCM.mutect2_snv.tsv',sep='\t',index_col=0)  # 467 samples have mutations
mutation = mutation.loc[mutation['filter']=='PASS',:]
burden = pd.read_csv('result/burden_stage3.txt',sep='\t',index_col=0).loc['burden',:].iloc[:-1]  # 472
burden.index = ['-'.join(sample.split('-')[0:4]) for sample in burden.index]
snaf.mutation_analysis(mode='compute',burden=burden,mutation=mutation,output='result/stage3_mutation.txt')


'''neoantigen analysis'''
snaf.downstream.analyze_neoantigens(freq_path='result/frequency_stage3_verbosity1_uid.txt',junction_path='result/burden_stage0.txt',total_samples=472,outdir='result',mers=None,fasta=False)
snaf.run_dash_T_antigen(input_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/result/shared_vs_unique_neoantigen_all.txt',
                        output_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/TCGA/SKCM/snaf_analysis/result')
snaf.downstream.plot_umap_neoantigen(df_path='result/mer9_umap_embed_df.txt',outdir='result')
        




























