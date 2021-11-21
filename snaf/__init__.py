from .snaf import snaf_configuration,NeoJunction,JunctionCountMatrixQuery
from .gtex import gtex_configuration
from .gtex_viewer import gtex_viewer_configuration
from .proteomics import *
from datetime import datetime,date

def initialize(exon_table,fasta,software_path,gtex_db,t_min=20,n_max=3):
    print('{} {} starting initialization'.format(date.today(),datetime.now().strftime('%H:%M:%S')))
    snaf_configuration(exon_table,fasta,software_path)
    gtex_configuration(gtex_db,t_min,n_max)
    gtex_viewer_configuration(gtex_db)
    print('{} {} finishing initialization'.format(date.today(),datetime.now().strftime('%H:%M:%S')))
