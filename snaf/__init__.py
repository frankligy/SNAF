from .snaf import snaf_configuration,NeoJunction,JunctionCountMatrixQuery,uid_to_coord, add_coord_frequency_table
from .gtex import gtex_configuration, tumor_specificity, add_tumor_specificity_frequency_table
from .gtex_viewer import gtex_viewer_configuration, gtex_visual_combine, gtex_visual_subplots
from .proteomics import *
from .downstream import *
from datetime import datetime,date
from .dash_app import run_dash_T_antigen,run_pweblogo
import os,sys

def initialize(db_dir,gtex_mode,software_path=None,binding_method=None,t_min=20,n_max=3,add_control=None):
    print('{} {} starting initialization'.format(date.today(),datetime.now().strftime('%H:%M:%S')))
    exon_table = os.path.join(db_dir,'Hs_Ensembl_exon_add_col.txt')
    fasta = os.path.join(db_dir,'Hs_gene-seq-2000_flank.fa')  # 2000 only affect the query_from_dict_fa function
    if gtex_mode == 'count':
        gtex_db = os.path.join(db_dir,'GTEx_junction_counts.h5ad')
    elif gtex_mode == 'psi':
        gtex_db = os.path.join(db_dir,'GTEx_junction_psi.h5ad')
    snaf_configuration(exon_table,fasta,software_path,binding_method)
    gtex_configuration(gtex_db,t_min,n_max,add_control)
    gtex_viewer_configuration(gtex_db)
    print('{} {} finishing initialization'.format(date.today(),datetime.now().strftime('%H:%M:%S')))

    # create a scratch folder
    scratch_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'scratch')
    if not os.path.exists(scratch_dir):
        os.mkdir(scratch_dir)



