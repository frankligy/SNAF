from .snaf import snaf_configuration,NeoJunction,JunctionCountMatrixQuery
from .gtex import gtex_configuration
from .gtex_viewer import gtex_viewer_configuration
from .binding import binding_configuration
from .proteomics import *
from datetime import datetime,date
from .dash_app import run_dash_app,run_pweblogo

def initialize(exon_table,fasta,gtex_db,software_path=None,binding_method=None,t_min=20,n_max=3):
    print('{} {} starting initialization'.format(date.today(),datetime.now().strftime('%H:%M:%S')))
    snaf_configuration(exon_table,fasta,software_path,binding_method)
    binding_configuration(binding_method)
    gtex_configuration(gtex_db,t_min,n_max)
    gtex_viewer_configuration(gtex_db)
    print('{} {} finishing initialization'.format(date.today(),datetime.now().strftime('%H:%M:%S')))



