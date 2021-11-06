from .snaf import snaf_configuration,NeoJunction
from .gtex import gtex_configuration

def initialize(exon_table,fasta,software_path,gtex_db):
    snaf_configuration(exon_table,fasta,software_path)
    gtex_configuration(gtex_db)
