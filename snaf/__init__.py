from .snaf import snaf_configuration,NeoJunction,JunctionCountMatrixQuery,uid_to_coord, add_coord_frequency_table, enhance_frequency_table
from .gtex import gtex_configuration, tumor_specificity, add_tumor_specificity_frequency_table
from .gtex_viewer import gtex_viewer_configuration, gtex_visual_combine, gtex_visual_subplots
from .proteomics import *
from .downstream import *
from datetime import datetime,date
from .dash_app import run_dash_T_antigen,run_pweblogo
import os,sys

def initialize(db_dir,gtex_mode,software_path=None,binding_method=None,t_min=20,n_max=3,add_control=None):
    '''
    Setting up global variable for running the program

    :param db_dir: string, the path to the database
    :param gtex_mode: string, either 'psi' or 'count'
    :param software_path: string or None, either the path to the netMHCpan4.1 executable or None (using MHCflurry)
    :param binding_method: string or None, either 'netMHCpan' or 'MHCflurry'
    :param t_min: int, the minimum number of read count the tumor sample should be larget than average number in normal database
    :param n_max: int, the maximum number of average read count normal database count have
    :param add_control: dataframe, the df containing the count matrix for additional matched control

    Example::

        # database directory (where you extract the reference tarball file)
        db_dir = '/user/ligk2e/download'
        # instantiate (if using netMHCpan)
        netMHCpan_path = '/user/ligk2e/netMHCpan-4.1/netMHCpan'
        snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='netMHCpan',software_path=netMHCpan_path)
        # instantiate (if not using netMHCpan)
        snaf.initialize(db_dir=db_dir,gtex_mode='count',binding_method='MHCflurry',software_path=None)

    '''
    print('{} {} starting initialization'.format(date.today(),datetime.now().strftime('%H:%M:%S')))
    exon_table = os.path.join(db_dir,'Hs_Ensembl_exon_add_col.txt')
    transcript_db = os.path.join(db_dir,'mRNA-ExonIDs.txt')
    fasta = os.path.join(db_dir,'Hs_gene-seq-2000_flank.fa')  # 2000 only affect the query_from_dict_fa function
    if gtex_mode == 'count':
        gtex_db = os.path.join(db_dir,'GTEx_junction_counts.h5ad')
    elif gtex_mode == 'psi':
        gtex_db = os.path.join(db_dir,'GTEx_junction_psi.h5ad')
    snaf_configuration(exon_table,transcript_db,db_dir,fasta,software_path,binding_method)
    gtex_configuration(gtex_db,t_min,n_max,add_control)
    gtex_viewer_configuration(gtex_db)
    print('{} {} finishing initialization'.format(date.today(),datetime.now().strftime('%H:%M:%S')))

    # create a scratch folder
    scratch_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'scratch')
    if not os.path.exists(scratch_dir):
        os.mkdir(scratch_dir)


def get_reduced_junction_matrix(pc,pea,frac=None,n=None):
    # preprocess the dataframe
    df = pd.read_csv(pc,sep='\t',index_col=0)
    df.index = [item.split('=')[0] for item in df.index]
    df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]
    print('original shape: {}'.format(df.shape))

    # filter to EventAnnotation file
    ee = [':'.join(item.split('|')[0].split(':')[1:]) for item in pd.read_csv(pea,sep='\t')['UID']]
    df = df.loc[df.index.isin(set(ee)),:]

    # for testing purpose
    if frac is not None:
        df = df.sample(frac=frac)
    if n is not None:
        df = df.sample(n=n)
    print('current shape: {}'.format(df.shape))
    return df




