from .snaf import snaf_configuration,NeoJunction,JunctionCountMatrixQuery
from .gtex import gtex_configuration

def initialize(exon_table,fasta,software_path,gtex_db,t_min=20,n_max=3):
    snaf_configuration(exon_table,fasta,software_path)
    gtex_configuration(gtex_db,t_min,n_max)
