#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import os
import sys
sys.path.insert(0,'/data/salomonis2/software')
import snaf
from snaf import surface
import pandas as pd
import numpy as np
import pickle
import subprocess
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from copy import deepcopy
import anndata as ad

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

# # get junction count matrix
# df = snaf.get_reduced_junction_matrix(pc='counts.original.txt',pea='Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt')

# # run SNAF
# netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
# db_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/data'
# tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
# gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
# add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}

# snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
# surface.initialize(db_dir=db_dir)

# # T antigen
# jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=30,add_control=add_control,outdir='result',filter_mode='maxmin')
# sample_to_hla = pd.read_csv('./sample_hla.txt',sep='\t',index_col=2)['HLA'].to_dict()
# hlas = [hla_string.split(',') for hla_string in df.columns.map(sample_to_hla)]
# jcmq.run(hlas=hlas,outdir='./result')
# snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')
# # proteome
# jcmq = snaf.JunctionCountMatrixQuery.deserialize('./result/after_prediction.p')
# for sample in jcmq.junction_count_matrix.columns:
#     print('polishing {} results'.format(sample))
#     jcmq.show_neoantigen_as_fasta(outdir='./fasta',name='neoantigen_{}.fasta'.format(sample),stage=2,verbosity=1,contain_uid=False,sample=sample)
#     snaf.proteomics.remove_redundant('./fasta/neoantigen_{}.fasta'.format(sample),'./fasta/neoantigen_{}_unique.fasta'.format(sample))
#     snaf.proteomics.compare_two_fasta(fa1_path='/data/salomonis2/LabFiles/Frank-Li/clinical/ovarian/MS/database/human_proteome_uniprot_9_10_mers_unique.fasta',
#                                       fa2_path='./fasta/neoantigen_{}_unique.fasta'.format(sample),outdir='./fasta',
#                                       write_unique2=True,prefix='{}_'.format(sample))

# # configure maxquant
# fasta_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ovarian/fasta'
# raw_dir = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ovarian/MS'
# root_dir = os.path.dirname(os.path.abspath(__file__))
# dic = {
#     'OvCa48':'SRR5933726',
#     'OvCa53':'SRR5933729',
#     'OvCa58':'SRR5933728',
#     'OvCa64':'SRR5933735',
#     'OvCa65':'SRR5933734',
#     'OvCa70':'SRR5933738',
#     'OvCa80':'SRR5947644',
#     'OvCa84':'SRR5947645',
#     'OvCa99':'SRR5947646',
#     'OvCa104':'SRR5947647',
#     'OvCa105':'SRR5933743',
#     'OvCa109':'SRR5933745',
#     'OvCa111':'SRR5933736',
#     'OvCa114':'SRR5933737'
# }

# for k,v in dic.items():
#     dbs = [os.path.join(fasta_dir,'{}_secondAligned.sortedByCoord.out.bed_unique2.fasta'.format(v))]
#     os.chdir(os.path.join(raw_dir,k))
#     inputs = subprocess.run("for file in *.raw; do echo $file; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
#     os.chdir(root_dir)
#     inputs = [os.path.join(raw_dir,k,inp) for inp in inputs]
#     outdir = os.path.join(raw_dir,k)
#     snaf.proteomics.set_maxquant_configuration(dbs=dbs,n_threads=20,inputs=inputs,enzymes=None,enzyme_mode=5,protein_fdr=1,peptide_fdr=0.05,site_fdr=1,outdir=outdir)



# specifically test candidates
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

control_bam_path_list = ['/data/salomonis-archive/BAMs/PublicDatasets/E-MTAB-2836-Grch38_Deep-Healthy-PanTissue/ERR315455_1.bam']
control_bai_path_list = ['/data/salomonis-archive/BAMs/PublicDatasets/E-MTAB-2836-Grch38_Deep-Healthy-PanTissue/ERR315455_1.bam.bai']

uid = 'ENSG00000196421:E2.1-E3.1'
sample = 'SRR5947644_secondAligned.sortedByCoord.out.bed'
region = 'chr20:64036326-64037927'
bam_path_list = ['/data/salomonis-archive/FASTQs/NCI-R01/SNAF_ovarian/SRR5947644_secondAligned.sortedByCoord.out.bam'] + control_bam_path_list
bai_path_list = ['/data/salomonis-archive/FASTQs/NCI-R01/SNAF_ovarian/SRR5947644_secondAligned.sortedByCoord.out.bam.bai'] + control_bai_path_list
sif_anno_path = '/data/salomonis2/software/ggsashimi'
outdir = 'Frank_inspection/sashimi'
bam_contig_rename = [False,True]
criterion=[('netMHCpan_el', 0, '<=', 2)]

def flank_chrom(chrom,offset):
    chrom = chrom.split('(')[0]
    sec1 = chrom.split(':')[0]
    sec2 = chrom.split(':')[1].split('-')[0]
    sec3 = chrom.split(':')[1].split('-')[1]
    new_sec2 = int(sec2) - offset[0]
    new_sec3 = int(sec3) + offset[1]
    assemble = '{}:{}-{}'.format(sec1,new_sec2,new_sec3)
    return assemble

# snaf.gtex_visual_combine_plotly(uid=uid,outdir='Frank_inspection',norm=False,tumor=df)
# snaf.gtex_visual_combine_plotly(uid=uid,outdir='Frank_inspection',norm=True,tumor=df)
# snaf.JunctionCountMatrixQuery.deserialize('result/after_prediction.p').visualize(uid=uid,sample=sample,outdir='Frank_inspection',criterion=criterion)
snaf.prepare_sashimi_plot(bam_path_list,bai_path_list,outdir,sif_anno_path,bam_contig_rename,query_region=flank_chrom(region,(2000,2000)),skip_copy=False, min_junction=1)
sys.exit('stop')

'''
NQDEDPLEV:
ENSG00000137434:E6.1-E7.1
chr6:10672636-10683187
SRR5947646 OvCa99
testis, ERR315351

KGPWYPLSL
ENSG00000196421:E2.1-E3.1
chr20:64036326-64037927
C20orf204
SRR5947644  OvCa80
brain, ERR315455
'''

# stack barplot
srr_to_id = pd.read_csv('meta.txt',sep='\t',index_col=0).squeeze().to_dict()
fig,ax = plt.subplots()
all_n = []
all_v = []
for i,(srr,id_) in enumerate(srr_to_id.items()):
    db = '{}_secondAligned.sortedByCoord.out.bed_unique2.fasta'.format(srr)
    n = subprocess.run(['wc','-l','./fasta/{}'.format(db)],stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[0].split(' ')[0]
    n = int(n)/2
    all_n.append(n)
    pep_path = 'MS/{}/combined/txt/peptides.txt'.format(id_)
    pep = pd.read_csv(pep_path,sep='\t',index_col=0)
    pep = pep.loc[pep['Proteins'].notna(),:]
    v = pep.shape[0]
    v = int(v)
    all_v.append(v)
    ax.bar(x=i,height=n,bottom=0,color='b')
    ax.bar(x=i,height=v,bottom=n-v,color='orange')
    ax.text(x=i,y=n+5,s=v,fontsize=5,ha='center')
ax.legend(handles=[Patch(color=i) for i in ['b','orange']],labels=['predicted','MS supported'],loc='upper left',bbox_to_anchor=(1,1),frameon=False)
ax.set_xticks(np.arange(len(srr_to_id)))
ax.set_xticklabels([item.replace('OvCa','P')for item in srr_to_id.values()],rotation=60,fontsize=8)
ax.set_xlabel('Ovarian Cancer Patients')
ax.set_ylabel('Number of Peptides')
plt.savefig('stack_barplot.pdf',bbox_inches='tight')
plt.close()

rate = np.array(all_v) / np.array(all_n)
'''
[77, 28, 18, 16, 12, 21, 160, 143, 43, 12, 25, 15, 29, 42]
[1140.0, 439.0, 739.0, 348.0, 579.0, 330.0, 988.0, 446.0, 732.0, 529.0, 843.0, 519.0, 502.0, 604.0]
[0.06754386 0.06378132 0.02435724 0.04597701 0.02072539 0.06363636
 0.16194332 0.3206278  0.05874317 0.02268431 0.02965599 0.02890173
 0.05776892 0.06953642]

 0.07399163277018352
'''
sys.exit('stop')







# # output supp1 table, MS results and other information
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS')
# all_patients = subprocess.run('for folder in *; do if [ -d $folder ]; then echo $folder; fi; done',shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
# os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis')
# patient_2_sample = pd.read_csv('sample_hla.txt',sep='\t',index_col=0)['sample_name'].to_dict()
# freq_df = pd.read_csv('result/frequency_stage2_verbosity1_uid_gene_coord_mean_mle.txt',sep='\t',index_col=0)
# aa_2_uid_first = {}
# uid_2_ts_mean = {}
# uid_2_ts_mle = {}
# for item,ts_mean,ts_mle in zip(*[freq_df.index,freq_df['tumor_specificity_mean'],freq_df['tumor_specificity_mle']]):
#     aa,uid = item.split(',')
#     aa_2_uid_first.setdefault(aa,[]).append(uid)
#     uid_2_ts_mean[uid] = ts_mean
#     uid_2_ts_mle[uid] = ts_mle
# aa_2_uid_second = {aa:','.join(uid) for aa,uid in aa_2_uid_first.items()}
# jcmq = snaf.JunctionCountMatrixQuery.deserialize('result/after_prediction.p')
# with pd.ExcelWriter('/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/supp1_table.xlsx') as writer:
#     for p in all_patients:
#         print(p)
#         sample = patient_2_sample[p]
#         df = pd.read_csv('result/frequency_stage2_verbosity1_uid.txt',sep='\t',index_col=0)
#         snaf.report_candidates(jcmq,df,sample=sample,outdir='result/candidates',criterion=[('netMHCpan_el',0,'<=',2),])
#         cand = pd.read_csv('result/candidates/T_antigen_candidates_{}.txt'.format(sample),sep='\t',index_col=0)
#         aa_2_binding = {}
#         aa_2_immuno = {}
#         for antigen,sub_df in cand.groupby(by='peptide'):
#             dict_binding = pd.Series(index=sub_df['hla'].values,data=sub_df['binding_affinity'].values).to_dict()
#             dict_immuno = pd.Series(index=sub_df['hla'].values,data=sub_df['immunogenicity'].values).to_dict()
#             aa_2_binding[antigen] = dict_binding
#             aa_2_immuno[antigen] = dict_immuno
#         os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/MS/{}/combined/txt'.format(p))
#         pep = pd.read_csv('peptides.txt',sep='\t',index_col=0)
#         col_uid = []
#         col_ts_mean = []
#         col_ts_mle = []
#         col_nmp = []
#         col_di = []
#         for aa,protein in zip(*[pep.index,pep['Proteins']]):
#             if pd.notna(protein):
#                 uid = aa_2_uid_second[aa]
#                 means = []
#                 mles = []
#                 for u in uid.split(','):
#                     means.append(str(uid_2_ts_mean[u]))
#                     mles.append(str(uid_2_ts_mle[u]))
#                 means = ','.join(means)
#                 mles = ','.join(mles)
#                 nmp = aa_2_binding[aa]
#                 di = aa_2_immuno[aa]
#             else:
#                 uid, means, mles, nmp, di = pd.NA, pd.NA, pd.NA, pd.NA, pd.NA
#             col_uid.append(uid)
#             col_ts_mean.append(means)
#             col_ts_mle.append(mles)
#             col_nmp.append(nmp)
#             col_di.append(di)
#         pep['uid'] = col_uid
#         pep['tumor_specificity_mean'] = col_ts_mean
#         pep['tumor_specificity_mle'] = col_ts_mle
#         pep['binding_netMHCpan'] = col_nmp
#         pep['immunogenicity_deepimmuno'] = col_di
#         pep.to_excel(writer,sheet_name=p)
#         os.chdir('/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis')

with pd.ExcelWriter('/data/salomonis2/LabFiles/Frank-Li/neoantigen/MS/schuster/RNA/snaf_analysis/supp1_table.xlsx', mode='a') as writer:
    # add look up table
    lookup = pd.read_csv('sample_hla.txt',sep='\t',index_col=0)
    lookup.to_excel(writer,sheet_name='lookup_table')
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
        'uid':'The parental NeoJunction of the amino acid sequence',
        'tumor_specificity_mean':'The average raw read counts for each parental NeoJunction',
        'tumor_specificity_mle':'The Tumor specificity score derived using Maximum Likelihood Estimation(MLE)',
        'binding_netMHCpan':'netMHCpan binding affinity of corresponding peptide-HLA complex',
        'immunogenicity_deepimmuno':'deepimmuno immunogenicity score of the corresponding peptide-HLA complex'
    }
    pd.Series(data=dic,name='description').to_excel(writer,sheet_name='description')




    














