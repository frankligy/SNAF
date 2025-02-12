#!/gpfs/data/yarmarkovichlab/Frank/test_snaf/test_snaf_env/bin/python3.7

import os,sys
import pandas as pd
import numpy as np
import snaf
import anndata as ad


# change the below two paths
db_dir = '/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/data'
netMHCpan_path = '/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/netMHCpan-4.1/netMHCpan' 

# initiate
tcga_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','tcga_matched_control_junction_count.h5ad'))
gtex_skin_ctrl_db = ad.read_h5ad(os.path.join(db_dir,'controls','gtex_skin_count.h5ad'))
add_control = {'tcga_control':tcga_ctrl_db,'gtex_skin':gtex_skin_ctrl_db}
df = pd.read_csv('counts.original.pruned.txt',sep='\t',index_col=0)
snaf.initialize(df=df,db_dir=db_dir,binding_method='netMHCpan',software_path=netMHCpan_path,add_control=add_control)
print('-------------pass initiate test----------------')

# test T antigen function
jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df,cores=8,add_control=add_control,outdir='result')
hlas = [['HLA-A*02:01','HLA-A*02:01','HLA-B*39:10','HLA-B*15:01','HLA-C*03:03','HLA-C*12:03'],
        ['HLA-A*02:01','HLA-A*01:01','HLA-B*40:01','HLA-B*52:01','HLA-C*03:04','HLA-C*12:02']]
jcmq.run(hlas=hlas,outdir='./result')
snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')
print('-------------pass T antigen function test----------------')

# test T antigen visual and server, for dash T viewer, please modify the full path accordingly
jcmq = snaf.JunctionCountMatrixQuery.deserialize('result/after_prediction.p')
uid = 'ENSG00000183856:E10.1-E12.7'
jcmq.visualize(uid=uid,sample='SRR5933726.Aligned.sortedByCoord.out.bed',outdir='./result')
dff = snaf.gtex_visual_combine(uid=uid,outdir='result',norm=False,tumor=df,group_by_tissue=False)
snaf.analyze_neoantigens(freq_path='result/frequency_stage2_verbosity1_uid.txt',junction_path='result/burden_stage0.txt',total_samples=2,outdir='result',mers=None,fasta=False)
print('-------------pass visualization test and please check your T antigen viewer following the prompt and URL, after checking, you can close this session, SNAF-T test done----------------')
snaf.run_dash_T_antigen(input_abs_path='/gpfs/data/yarmarkovichlab/Frank/test_snaf/result/shared_vs_unique_neoantigen_all.txt')


# test B antigen, you can comment out "test T antigen viusual and server" section
# change the software_path accordingly
# change the python_executable accordingly
from snaf import surface
surface.initialize(db_dir=db_dir)
print('-------------pass B antigen initiate test----------------')

df = pd.read_csv('counts.original.pruned.txt',sep='\t',index_col=0)
membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df,add_control=add_control,outdir='result/surface')
# short_read mode
surface.run(uids=membrane_tuples,outdir='result/surface',prediction_mode='short_read',
            gtf=None,
            tmhmm=True,software_path='/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/tmhmm-2.0c/bin/tmhmm') 
surface.generate_full_results(outdir='result/surface',mode='short_read',
                              freq_path='result/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt',
                              validation_gtf=os.path.join(db_dir,'2021UHRRIsoSeq_SQANTI3_filtered.gtf'))
print('-------------pass B antigen short read mode test----------------')

# long_read mode
surface.run(uids=membrane_tuples,outdir='result/surface',prediction_mode='long_read',
            gtf=os.path.join(db_dir,'2021UHRRIsoSeq_SQANTI3_filtered.gtf'),
            tmhmm=True,software_path='/gpfs/data/yarmarkovichlab/Frank/SNAF_ecosystem/tmhmm-2.0c/bin/tmhmm')
surface.generate_full_results(outdir='result/surface',mode='long_read',
                              freq_path='result/frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt',
                              validation_gtf=None)
print('-------------pass B antigen short long mode test----------------')
print('-------------SNAF-B viewer will be launched, follow prompt, once finished, you can close the session----------------')

surface.run_dash_B_antigen(pkl='result/surface/surface_antigen_lr.p',candidates='result/surface/candidates_3_lr_None.txt',prediction_mode='long_read',
                           python_executable='/gpfs/data/yarmarkovichlab/Frank/test_snaf/test_snaf_env/bin/python3.7')


