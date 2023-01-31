from .snaf import snaf_configuration,NeoJunction,JunctionCountMatrixQuery,uid_to_coord, add_coord_frequency_table, enhance_frequency_table
from .gtex import gtex_configuration, tumor_specificity, add_tumor_specificity_frequency_table
from .gtex_viewer import gtex_viewer_configuration, gtex_visual_combine, gtex_visual_subplots, gtex_visual_combine_plotly, gtex_visual_per_tissue_count
from .proteomics import *
from .downstream import *
from datetime import datetime,date
from .dash_app import run_dash_T_antigen,run_pweblogo
import os,sys

def initialize(df,db_dir,gtex_mode='count',software_path=None,binding_method=None,t_min=20,n_max=3,
               normal_cutoff=5, tumor_cutoff=20, normal_prevalance_cutoff=0.01, tumor_prevalance_cutoff=0.1, add_control=None):
    '''
    Setting up global variable for running the program

    :param df: dataframe, the currently tested cancer cohort dataframe, this will limit number of junctions being loaded into memory for control database
    :param db_dir: string, the path to the database
    :param gtex_mode: string, either 'psi' or 'count'
    :param software_path: string or None, either the path to the netMHCpan4.1 executable or None (using MHCflurry)
    :param binding_method: string or None, either 'netMHCpan' or 'MHCflurry'
    :param t_min: int, the minimum number of read count the tumor sample should be larget than average number in normal database
    :param n_max: int, the maximum number of average read count normal database count have
    :param normal_cutoff: int, below which read count we consider a junction is not expressed in normal tissue
    :param tumor_cutoff: int, above which read count we consider a junction is expressed in tumor tissue
    :param normal_prevalance_cutoff: float, if below this fraction, we consider a junction is not present in normal tissue at an appreciable amount
    :param tumor_prevalance_cutoff: float, if above this fraction, we consider a junction is present in tumor tissue at an appreciable amount
    :param add_control: None or a dictionary containing additional controls, additional controls can a dataframe or anndata, for instance, if adding two controls asides from internal GTEx,
                        {'tcga_matched_control':df,'gtex_skin':adata}, when added database contains multiple tissue types, it is recommended to pass as a AnnData with the tissue
                        information stored as adata.var['tissue'], otherwise, if passed as a dataframe, or samples will be recognized as a single tissue type, it will affect tumor specifcity
                        score calculation if using MLE or bayesian hierarchical models.

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
    exon_table = os.path.join(db_dir,'Alt91_db','Hs_Ensembl_exon_add_col.txt')
    transcript_db = os.path.join(db_dir,'Alt91_db','mRNA-ExonIDs.txt')
    fasta = os.path.join(db_dir,'Alt91_db','Hs_gene-seq-2000_flank.fa')  # 2000 only affect the query_from_dict_fa function
    if gtex_mode == 'count':
        gtex_db = os.path.join(db_dir,'controls','GTEx_junction_counts.h5ad')
        # gtex_db = '/data/salomonis2/LabFiles/Frank-Li/neoantigen/revision/ts/gene/gtex_count.h5ad'
    elif gtex_mode == 'psi':
        gtex_db = os.path.join(db_dir,'controls','GTEx_junction_psi.h5ad')
    snaf_configuration(exon_table,transcript_db,db_dir,fasta,software_path,binding_method)
    adata = gtex_configuration(df,gtex_db,t_min,n_max,normal_cutoff, tumor_cutoff, normal_prevalance_cutoff, tumor_prevalance_cutoff, add_control)
    gtex_viewer_configuration(adata)
    print('{} {} finishing initialization'.format(date.today(),datetime.now().strftime('%H:%M:%S')))

    # create a scratch folder
    scratch_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),'scratch')
    if not os.path.exists(scratch_dir):
        os.mkdir(scratch_dir)


def get_reduced_junction_matrix(pc,pea,frac=None,n=None,samples=None):
    # preprocess the dataframe
    df = pd.read_csv(pc,sep='\t',index_col=0)
    df.index = [item.split('=')[0] for item in df.index]
    df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]
    print('original shape: {}'.format(df.shape))
    if samples is not None:
        print('{} samples have HLA type info, will subset to these samples'.format(len(samples)))

    # filter to EventAnnotation file
    ee = [':'.join(item.split('|')[0].split(':')[1:]) for item in pd.read_csv(pea,sep='\t')['UID']]
    df = df.loc[df.index.isin(set(ee)),:]

    # for testing purpose
    if frac is not None:
        df = df.sample(frac=frac)
    if n is not None:
        df = df.sample(n=n)

    # if need filtering out or rearranging samples
    if samples is not None:
        df = df.loc[:,df.columns.isin(samples)] 
    print('current shape: {}'.format(df.shape))
    return df


def remove_trailing_coord(count,sep='\t'):
    df = pd.read_csv(count,sep=sep,index_col=0)
    df.index = [item.split('=')[0] for item in df.index]
    df = df.loc[np.logical_not(df.index.duplicated()).tolist(),:]
    return df




