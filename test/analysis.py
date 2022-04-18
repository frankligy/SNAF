#!/data/salomonis2/LabFiles/Frank-Li/neoantigen/step2/snaf_test_env/bin/python3.7

# T antigen
import os,sys
import pandas as pd
import numpy as np
import snaf


db_dir = '/data/salomonis2/LabFiles/Frank-Li/refactor/download'
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='netMHCpan',software_path=netMHCpan_path)
df = pd.read_csv('counts.original.pruned.txt',sep='\t',index_col=0)
jcmq = snaf.JunctionCountMatrixQuery(junction_count_matrix=df)
hlas = [['HLA-A*02:01','HLA-A*02:01','HLA-B*39:10','HLA-B*15:01','HLA-C*03:03','HLA-C*12:03'],
        ['HLA-A*02:01','HLA-A*01:01','HLA-B*40:01','HLA-B*52:01','HLA-C*03:04','HLA-C*12:02']]
jcmq.run(hlas=hlas,outdir='./result')
snaf.JunctionCountMatrixQuery.generate_results(path='./result/after_prediction.p',outdir='./result')

jcmq = snaf.JunctionCountMatrixQuery.deserialize('result/after_prediction.p')
jcmq.visualize(uid='ENSG00000137868:E7.2-E8.1',sample='SRR5933726.Aligned.sortedByCoord.out.bed',outdir='./result')
snaf.analyze_neoantigens(freq_path='result/frequency_stage2_verbosity1_uid.txt',junction_path='result/burden_stage0.txt',total_samples=2,outdir='result',mers=None,fasta=False)
snaf.run_dash_T_antigen(input_abs_path='/data/salomonis2/LabFiles/Frank-Li/neoantigen/step2/result/shared_vs_unique_neoantigen_all.txt')


# B antigen
import snaf
import pandas as pd
db_dir = '/data/salomonis2/LabFiles/Frank-Li/refactor/download'
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='netMHCpan',software_path=netMHCpan_path)
from snaf import surface
surface.initialize(db_dir=db_dir)

df = pd.read_csv('counts.original.pruned.txt',sep='\t',index_col=0)
membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df)
surface.run(membrane_tuples,outdir='result',tmhmm=True,software_path='/data/salomonis2/LabFiles/Frank-Li/python3/TMHMM/tmhmm-2.0c/bin/tmhmm')
surface.generate_results(pickle_path='./result/surface_antigen.p',outdir='result',strigency=3,gtf=None)
surface.run_dash_B_antigen(pkl='result/surface_antigen.p',candidates='result/candidates.txt',
                           python_executable='/data/salomonis2/LabFiles/Frank-Li/neoantigen/step2/snaf_test_env/bin/python3.7')


# tumor specificity
import snaf
import pandas as pd
db_dir = '/data/salomonis2/LabFiles/Frank-Li/refactor/download'
netMHCpan_path = '/data/salomonis2/LabFiles/Frank-Li/refactor/external/netMHCpan-4.1/netMHCpan'
snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='netMHCpan',software_path=netMHCpan_path)

df = pd.read_csv('counts.original.pruned.txt',sep='\t',index_col=0)
snaf.gtex_visual_combine('ENSG00000001561:E3.1-E4.1',norm=True,outdir='result',tumor=df)
snaf.gtex_visual_subplots('ENSG00000001561:E3.1-E4.1',norm=True,outdir='result')
df = pd.read_csv('result/frequency_stage2_verbosity1_uid.txt',sep='\t',index_col=0)
df = snaf.add_gene_symbol_frequency_table(df=df,remove_quote=True)
print(df)
df = snaf.add_coord_frequency_table(df=df,remove_quote=False)
print(df)