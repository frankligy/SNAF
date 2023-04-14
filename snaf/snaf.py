#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
from functools import reduce
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib as mpl
import anndata as ad
import seaborn as sns
from scipy import stats
from scipy.optimize import minimize
import re
import requests
import xmltodict
import multiprocessing as mp
import functools
from tqdm import tqdm
from itertools import compress
from ast import literal_eval
from bisect import bisect, bisect_left


# for biopython, pip install biopython
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio.Seq import Seq

# other modules
from .binding import run_netMHCpan,run_MHCflurry
from .deepimmuno import run_deepimmuno
from .data_io import *
from .visualize import *
from .gtex import *
from .gtex_viewer import *
from .downstream import *


import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


# configuration
def snaf_configuration(exon_table,transcript_db,db_dir,fasta,software_path_arg=None,binding_method_arg=None):
    global dict_exonCoords
    global dict_fa
    global software_path
    global binding_method
    global dict_exonlist
    global dict_start_codon
    global phase_inferer_gtf_dict
    dict_exonCoords = exonCoords_to_dict(exon_table)
    dict_exonlist = construct_dict_exonlist(transcript_db)
    dict_fa = fasta_to_dict(fasta)
    software_path = software_path_arg
    binding_method = binding_method_arg

    df_start_codon = pd.read_csv(os.path.join(db_dir,'Alt91_db','df_start_codon.txt'),sep='\t',index_col=0)
    df_start_codon['start_codon'] = [literal_eval(item) for item in df_start_codon['start_codon']]
    df_start_codon['non_redundant'] = [literal_eval(item) for item in df_start_codon['non_redundant']]
    dict_start_codon = df_start_codon['start_codon'].to_dict()

    phase_inferer_gtf_dict = process_gtf(os.path.join(db_dir,'Alt91_db','Homo_sapiens.GRCh38.91.gtf'))

    


def process_gtf(gtf):
    # any official ensembl release gtf format
    '''
    gtf_dict['ENSG00000186092']

        {'ENST00000641515': [(65419, 65433), (65520, 65573), (69037, 71585)],
        'ENST00000335137': [(69055, 70108)]}
    '''
    gtf_dict = {}
    with open(gtf,'r') as f:
        for line in f:
            try:
                chrom, source, typ, start, end, score, strand, phase, attrs = line.rstrip('\n').split('\t')
            except ValueError:  # the header line
                continue
            if typ == 'gene':
                ensg = attrs.split(';')[0].split(' ')[1].strip('"')
                gtf_dict[ensg] = {}
            elif typ == 'transcript':
                enst = attrs.split(';')[2].split(' ')[2].strip('"')
                gtf_dict[ensg][enst] = []
            elif typ == 'exon':
                gtf_dict[ensg][enst].append((int(start),int(end)))
    return gtf_dict

def get_support_phase(ensg, coord_first_exon_last_base, pssc, strand, length_first):
    coord_first_exon_last_base, pssc, length_first = int(coord_first_exon_last_base), int(pssc), int(length_first)
    all_trans = phase_inferer_gtf_dict[ensg]
    supports = []
    if strand == '+':
        for enst, tran in all_trans.items():
            lis = []
            [lis.extend(list(exon)) for exon in tran]
            pssc_insert_pos = bisect(lis, pssc)
            coord_first_exon_first_base = coord_first_exon_last_base - length_first + 1
            junction_insert_pos = bisect(lis, coord_first_exon_first_base)

            if junction_insert_pos % 2 == 1 and pssc_insert_pos % 2 == 1:
                if junction_insert_pos == pssc_insert_pos:
                    n_bases = coord_first_exon_first_base - pssc + 1
                elif junction_insert_pos > pssc_insert_pos:
                    start_exon_index = (pssc_insert_pos - 1) // 2
                    end_exon_index = (junction_insert_pos - 1) // 2
                    n_bases = 0
                    for i, exon in enumerate(tran):
                        if i == start_exon_index:
                            n_bases += (exon[1] - pssc + 1)
                        elif i == end_exon_index:
                            n_bases += (coord_first_exon_first_base - exon[0] + 1)
                        else:
                            if i > start_exon_index and i < end_exon_index:
                                n_bases += (exon[1] - exon[0] + 1)
                else:
                    continue
            else:
                continue


            remainder = n_bases % 3
            '''
            *   *   *   *   |*   *   *    *    *
                            |6   7   8    9    10

            remainder means a fragment including the first base in the first exon, how many left. 1 left means 6 should be the first base of new codon
            2 left means 6 should be 6 is the second in last codon, 7 will be the last base in the last codon, have to start with 8
            0 left means 6 should be the last base in the last codon, so that we start with 7
            '''
            if remainder == 1:
                phase = 0
            elif remainder == 2:
                phase = 2
            elif remainder == 0:
                phase = 1
            supports.append((phase, pssc, enst,strand))

    elif strand == '-':
        for enst, tran in all_trans.items():
            tran.sort(
                key=lambda x: x[0])  # because negative strand looks like [(5,7),(2,4)], we change it to [(2,4),(5,7)]
            lis = []
            [lis.extend(list(exon)) for exon in tran]
            pssc_insert_pos = bisect_left(lis, pssc)
            coord_first_exon_first_base = coord_first_exon_last_base + length_first - 1
            junction_insert_pos = bisect_left(lis, coord_first_exon_first_base)

            if junction_insert_pos % 2 == 1 and pssc_insert_pos % 2 == 1:
                if junction_insert_pos == pssc_insert_pos:
                    n_bases = pssc - coord_first_exon_first_base + 1
                elif junction_insert_pos < pssc_insert_pos:
                    start_exon_index = (pssc_insert_pos - 1) // 2
                    end_exon_index = (junction_insert_pos - 1) // 2
                    n_bases = 0
                    for i, exon in enumerate(tran):
                        if i == start_exon_index:
                            n_bases += (pssc - exon[0] + 1)
                        elif i == end_exon_index:
                            n_bases += (exon[1] - coord_first_exon_first_base + 1)
                        else:
                            if i > end_exon_index and i < start_exon_index:
                                n_bases += (exon[1] - exon[0] + 1)
                else:
                    continue
            else:
                continue



            remainder = n_bases % 3  # same idea as above forward strand, see explanation above
            if remainder == 1:
                phase = 0
            elif remainder == 2:
                phase = 2
            elif remainder == 0:
                phase = 1
            supports.append((phase, pssc, enst, strand))

    return supports   



    
def is_in_db(valid):
    # provide a list with uid without repeats
    mapping = {}
    for uid in tqdm(valid):
        ensg = uid.split(':')[0]
        exons = ':'.join(uid.split(':')[1:])
        if '_' in exons or 'U' in exons or 'ENSG' in exons or 'I' in exons:
            mapping[uid] = False
        else:
            exonlist = dict_exonlist[ensg]
            exonstring = '|'.join(exonlist)
            e1,e2 = exons.split('-')
            pattern1 = re.compile(r'^{}\|{}\|'.format(e1,e2))  # ^E1.1|E2.3|
            pattern2 = re.compile(r'\|{}\|{}$'.format(e1,e2))  # |E1.1|E2.3$
            pattern3 = re.compile(r'\|{}\|{}\|'.format(e1,e2)) # |E1.1|E2.3|
            if re.search(pattern3,exonstring) or re.search(pattern2,exonstring) or re.search(pattern1,exonstring):   # as long as match one pattern, should be eliminated
                mapping[uid] = True
            else:
                mapping[uid] = False
    return mapping








class JunctionCountMatrixQuery():

    '''
    Instantiate the JunctionCountMatrixQuery class

    :param junction_count_matrix: The pandas dataframe for the junction_count_matrix, outputted by AltAnalyze
    :param cores: int, how many cores you'd like to use, if None, then will automatically detect the maximum amount of cores in the OS
    :param add_control: None or a dictionary containing additional controls, additional controls can a dataframe or anndata, for instance, if adding two controls asides from internal GTEx,
                        {'tcga_matched_control':df,'gtex_skin':adata}, when added database contains multiple tissue types, it is recommended to pass as a AnnData with the tissue
                        information stored as adata.var['tissue'], otherwise, if passed as a dataframe, or samples will be recognized as a single tissue type, it will affect tumor specifcity
                        score calculation if using MLE or bayesian hierarchical models.
    :param not_in_db: boolean, whether to remove junctions that are present in any Ensembl documented transcript, remember some of the documented transcript in
                      Ensembl are tumor-specific as well, doing so may remove some bona fide hits. But good for reducing number of neoantigens for validation.
    :param outdir: string, the output folder for storing results
    :param filter_mode: string, either 'maxmin' or 'prevalance', this correspond to the parameters that will be used in snaf.initialize, if maxmin, use t_min and n_max, 
                        if prevalance, use normal and tumor cutoff

    Example::

        jcmq = JunctionCountMatrixQuery(junction_count_matrix=df,cores=50,add_control=control_df,not_in_db=True)
        
    '''

    def __init__(self,junction_count_matrix,cores=None,add_control=None,not_in_db=False,outdir='.',filter_mode='maxmin'):
        self.junction_count_matrix = junction_count_matrix
        if cores is None:
            cores = mp.cpu_count()
        self.cores = cores
        if not_in_db:
            self.get_neojunctions(add_control=add_control,dict_exonlist=dict_exonlist,outdir=outdir,filter_mode=filter_mode)
        else:
            self.get_neojunctions(add_control=add_control,dict_exonlist=None,outdir=outdir,filter_mode=filter_mode)
        

    def __str__(self):
        try:
            len_translated = len(self.translated)
        except AttributeError:
            len_translated = None
        try:
            shape_cond_subset_df = self.cond_subset_df.shape
        except AttributeError:
            shape_cond_subset_df = None
        try:
            len_results = len(self.results)
        except AttributeError:
            len_results = None
        return 'junction_count_matrix: {}\n'\
               'cores: {}\n'\
               'valid: {}\n'\
               'invalid: {}\n'\
               'cond_df: {}\n'\
               'subset: {}\n'\
               'translated: list of {} nj objects\n'\
               'cond_subset_df: {}\n'\
               'results: list of length {}'.format(self.junction_count_matrix.shape,self.cores,len(self.valid),len(self.invalid),
                                               self.cond_df.shape,self.subset.shape,len_translated,shape_cond_subset_df,len_results)
    
    def get_neojunctions(self,add_control,dict_exonlist,outdir,filter_mode):
        self.valid, self.invalid, self.cond_df = multiple_crude_sifting(self.junction_count_matrix,add_control,dict_exonlist,outdir,filter_mode)
        self.subset = self.junction_count_matrix.loc[self.valid,:]
        self.cond_subset_df = self.cond_df.loc[self.valid,:]

    def get_neojunction_info(self,event):  # more for B cell antigen, since we design it in a way that B cell pipeline completely rely on T pipeline for GTEx check
        ed = self.subset.loc[event,:].to_dict()
        freq = np.count_nonzero(np.array(list(ed.values())))/len(ed)
        return ed,freq


    @staticmethod
    def get_membrane_tuples(df,**kwargs):
        '''
        this function is used by SurfaceAntigen pipeline to filter out splicing evnets that are not tumor-specific and also compute
        useful informations for each membrane protein.

        :param df: pandas dataframe, the junction count matrix
        :param **kwargs: will be passed to __init__ of JunctionCountMatrixQuery class function

        :return membrane_tuples: a list, in which each item is a tuple (uid,mean_gtex,df_gtex,ed,freq). 

            * uid: the uid of the splicing evnet
            * mean_gtex: the mean raw read count across GTEx
            * df_gtex: a dataframe of gtex sample name, tissue types and read count values
            * ed: expression dictionary for tumor samples
            * freq: expression frequency for tumor samples

        Example::

            membrane_tuples = snaf.JunctionCountMatrixQuery.get_membrane_tuples(df)

        '''
        from .surface import filter_to_membrane_protein
        jcmq = JunctionCountMatrixQuery(junction_count_matrix=df,**kwargs)
        print(jcmq)
        neojunctions = jcmq.valid
        membrane_uid = filter_to_membrane_protein(neojunctions)
        membrane_tuples = []
        for uid in membrane_uid:
            mean_gtex, df_gtex = tumor_specificity(uid,method='mean',return_df=True)
            ed, freq = jcmq.get_neojunction_info(uid)
            membrane_tuples.append((uid,mean_gtex,df_gtex,ed,freq))
        return membrane_tuples



    @staticmethod    
    def split_df_to_chunks(df,cores=None):
        df_index = np.arange(df.shape[0])
        if cores is None:
            cores = mp.cpu_count()
        sub_indices = np.array_split(df_index,cores)
        sub_dfs = [df.iloc[sub_index,:] for sub_index in sub_indices]
        return sub_dfs

    @staticmethod
    def split_array_to_chunks(array,cores=None):
        if not isinstance(array,list):
            raise Exception('split_array_to_chunks function works for list, not ndarray')
        array_index = np.arange(len(array))
        if cores is None:
            cores = mp.cpu_count()
        sub_indices = np.array_split(array_index,cores)
        sub_arrays = []
        for sub_index in sub_indices:
            item_in_group = []
            for i in sub_index:
                item_in_group.append(array[i])
            sub_arrays.append(item_in_group)
        return sub_arrays
        
            
    @staticmethod
    def each_chunk_func(input_,kind,hlas=None,strict=False,sub_cond_df=None,binding_method=None):
        if kind == 1:
            nj_list = []
            for uid in input_.index:
                nj = NeoJunction(uid=uid,count=50,check_gtex=False)
                nj.detect_type()
                nj.retrieve_junction_seq()
                nj.in_silico_translation(strict=strict)    
                nj_list.append(nj) 
        elif kind == 2:   # out of development
            nj_list = []
            for uid in input_.index:
                nj = NeoJunction(uid=uid,count=50,check_gtex=False)
                nj.detect_type()
                nj.retrieve_junction_seq()
                nj.in_silico_translation()    
                nj.binding_prediction()
                nj.immunogenicity_prediction()
                nj_list.append(nj)  
        elif kind == 3:  # currently recommended
            nj_list = []
            assert len(input_) == sub_cond_df.shape[0]
            for nj,index in tqdm(zip(input_,range(sub_cond_df.shape[0])),total=len(input_)):
                cond_row = sub_cond_df.iloc[index].tolist()
                hlas_row = list(compress(hlas,cond_row))
                combined_unique_hlas = reduce(lambda a,b:list(set(a+b)),hlas_row)
                try:
                    nj.binding_prediction(hlas=combined_unique_hlas,binding_method=binding_method)
                except Exception as e:
                    if str(e) == 'Already no candidates after in-silico translation':
                        nj_list.append(None)
                        continue
                    else:
                        raise Exception('binding prediction error: {}'.format(e))
                try:
                    nj.immunogenicity_prediction()
                except Exception as e:
                    if str(e) == 'Already no candidates after binding prediction':
                        nj_list.append(None)
                        continue               
                    else:
                        raise Exception('immunogenic prediction error: {}'.format(e))     
                nj_list.append(nj) 
        elif kind == 4:
            nj_list = []
            for nj in input_:
                if nj is None:
                    nj_list.append(None)
                    continue
                try:
                    nj.binding_prediction(hlas=hlas)
                except Exception as e:
                    nj_list.append(None)
                    continue
                try:
                    nj.immunogenicity_prediction()
                except Exception as e:
                    nj_list.append(None)
                    continue
                nj_list.append(nj)                 
        return nj_list

    def run(self,hlas,strict=False,outdir='.',name='after_prediction.p'):
        '''
        main function to run mhc bound peptide T antigen pipeline

        :param hlas: list, each item is a list of 6 HLA types associated with each sample, the order of this list must be consistent with the sample name in the df column
        :param strict: boolean, default is False, determine whether to filter out peptide which doesn't have canonical start codon support
        :param outdir: string, the path where all the results will go into
        :param name: string, the name of the generated pickle result object

        Example::

            jcmq.run(hlas=hlas,outdir='result',name='after_prediction.p')
        '''
        self.parallelize_run(kind=1,strict=strict)
        print(self)
        self.parallelize_run(kind=3,hlas=hlas)
        self.serialize(outdir=outdir,name=name)

    @staticmethod
    def generate_results(path,outdir,criterion=None):
        '''
        wrapper function to automatically generate all necessary output

        :param path: string, where the pickle result file lie on the disk
        :param outdir: string, where to generate all the necessary result

        Example::

            JunctionCountMatrixQuery.generate_results(path='result/after_prediction.p',outdir='result')

        '''

        jcmq = JunctionCountMatrixQuery.deserialize(name=path)
        stage0_compatible_results(jcmq,outdir=outdir)
        for stage in [3,2,1]:
            jcmq.show_neoantigen_burden(outdir=outdir,name='burden_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False,criterion=criterion)
            jcmq.show_neoantigen_frequency(outdir=outdir,name='frequency_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False,plot=True,plot_name='frequency_stage{}.pdf'.format(stage),criterion=criterion)
            jcmq.show_neoantigen_frequency(outdir=outdir,name='frequency_stage{}_verbosity1_uid.txt'.format(stage),stage=stage,verbosity=1,contain_uid=True,plot=False,criterion=criterion)
            # add additional attributes
            df = pd.read_csv(os.path.join(outdir,'frequency_stage{}_verbosity1_uid.txt'.format(stage)),sep='\t',index_col=0)
            enhance_frequency_table(df,True,True,outdir,'frequency_stage{}_verbosity1_uid_gene_symbol_coord_mean_mle.txt'.format(stage))
            # report candidates
            if stage == 3:
                dff = pd.read_csv(os.path.join(outdir,'frequency_stage{}_verbosity1_uid_gene_symbol_coord_mean_mle.txt'.format(stage)),sep='\t',index_col=0)
                for sample in tqdm(jcmq.junction_count_matrix.columns,total=jcmq.junction_count_matrix.shape[1]):
                    report_candidates(jcmq,dff,sample,os.path.join(outdir,'T_candidates'),True)
                print('concatenating all T antigen files into one & indicate whether in AltAnalyze database')
                df_list = []
                for sample in jcmq.junction_count_matrix.columns:
                    df_list.append(pd.read_csv(os.path.join(outdir,'T_candidates','T_antigen_candidates_{}.txt'.format(sample)),sep='\t',index_col=0))
                final_df = pd.concat(df_list,axis=0)
                all_uid = final_df['uid'].unique().tolist()
                mapping = is_in_db(all_uid)
                final_df['in_db'] = final_df['uid'].map(mapping).values
                final_df.to_csv(os.path.join(outdir,'T_candidates','T_antigen_candidates_all.txt'),sep='\t')

        # add additional attributes to stage0
        df = pd.read_csv(os.path.join(outdir,'frequency_stage0.txt'),sep='\t',index_col=0)
        df.index = [','.join([item,item]) for item in df.index]
        df.to_csv(os.path.join(outdir,'frequency_stage0_verbosity1_uid.txt'),sep='\t')
        df = pd.read_csv(os.path.join(outdir,'frequency_stage0_verbosity1_uid.txt'),sep='\t',index_col=0)
        enhance_frequency_table(df,True,True,outdir,'frequency_stage0_verbosity1_uid_gene_symbol_coord_mean_mle.txt')        

    def parallelize_run(self,kind,hlas=None,strict=False):
        pool = mp.Pool(processes=self.cores)
        if kind == 1 or kind == 2:
            sub_dfs = JunctionCountMatrixQuery.split_df_to_chunks(self.subset,self.cores)
            r = [pool.apply_async(func=JunctionCountMatrixQuery.each_chunk_func,args=(sub_df,kind,None,strict,)) for sub_df in sub_dfs]
            pool.close()
            pool.join()
            results = []
            for collect in r:
                result = collect.get()
                results.extend(result)
            self.translated = results
        elif kind == 3:
            sub_arrays = JunctionCountMatrixQuery.split_array_to_chunks(self.translated,self.cores)
            sub_cond_dfs = JunctionCountMatrixQuery.split_df_to_chunks(self.cond_subset_df,self.cores)
            r = [pool.apply_async(func=JunctionCountMatrixQuery.each_chunk_func,args=(sub_array,kind,hlas,False,sub_cond_df,binding_method,)) for sub_array,sub_cond_df in zip(sub_arrays,sub_cond_dfs)]
            pool.close()
            pool.join()
            results = []
            for collect in r:
                result = collect.get()
                results.extend(result)
            self.results = (results,hlas)


        elif kind == 4: # heterogeneous hlas for each sample
            self.cond_subset_df = self.cond_df.loc[self.valid,:]
            need_to_predict_across_samples = []  # [[nj,None,nj,nj,None...],]
            for i,(label,content) in enumerate(self.cond_subset_df.iteritems()):
                need_to_predict = []
                number_need = 0
                for item in zip(self.translated,content.values):
                    if item[1]:
                        need_to_predict.append(item[0])
                        number_need += 1
                    else:
                        need_to_predict.append(None)
                need_to_predict_across_samples.append(need_to_predict)
                print('{} has {} junctions to proceed predictions, HLA are {}'.format(label,number_need,hlas[i]))
            print('Binding and immunogenicity prediction starts')
            r = [pool.apply_async(func=JunctionCountMatrixQuery.each_chunk_func,args=(need_to_predict,kind,hlas,)) for need_to_predict,hlas in zip(need_to_predict_across_samples,hlas)]
            pool.close()
            pool.join()
            results = []
            for collect in r:
                result = collect.get()
                results.append(result)
            self.results = results 

    def serialize(self,outdir,name):
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        with open(os.path.join(outdir,name),'wb') as f:
            pickle.dump(self,f,protocol=pickle.HIGHEST_PROTOCOL)  

    @staticmethod
    def deserialize(name):
        with open(name,'rb') as f:
            jcmq = pickle.load(f)
        return jcmq

    def visualize(self,uid,sample,outdir,tumor=False,criterion=[('netMHCpan_el',0,'<=',2),('deepimmuno_immunogenicity',1,'==','True'),]):
        '''
        Visualize certain Neojunction in certain sample

        :param uid: string, the uid for the event you want to check
        :param sample: string, the sample name
        :param outdir: string, where to deposite the figure
        :param tumor: bool, whether to show the expression level in tumor sample as well, default is not
        :param criterion: a nested tuple, if none of the neoantigens are bound and immunogenic, no candidates will be plotted, so you can relax the criterion a bit based on the format specified.

        Example::

            jcmq.visualize(uid='ENSG00000167291:E38.6-E39.1',sample='TCGA-DA-A1I1-06A-12R-A18U-07.bed',outdir='./result')
        '''
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if sample is not None:
            row_index = self.subset.index.tolist().index(uid)
            col_index = self.subset.columns.tolist().index(sample)
            results = self.results[0]
            hlas = self.results[1]
            nj = deepcopy(results[row_index])
            nj.enhanced_peptides = nj.enhanced_peptides.filter_based_on_hla(selected_hla=hlas[col_index])
            nj.visualize(outdir,'{}_{}.pdf'.format(uid.replace(':','_'),sample),criterion=criterion)
        if tumor:
            # in tumor sample
            fig,ax = plt.subplots()
            counts = self.junction_count_matrix.loc[uid,:]
            ax.plot(np.arange(len(counts)),counts.values,marker='o',markersize=3,color='k',markerfacecolor='r',markeredgecolor='r',linestyle='--')
            ax.set_xticks(np.arange(len(counts)))
            ax.set_xticklabels(counts.index.values,fontsize=1,rotation=90)
            ax.set_title('{}_tumor'.format(uid),fontsize=8)
            ax.set_ylim(bottom=-0.05)
            ax.set_ylabel('counts')
            plt.savefig(os.path.join(outdir,'{}_tumor.pdf'.format(uid.replace(':','_'))))
            plt.close()

    def show_neoantigen_burden(self,outdir,name,stage,verbosity,contain_uid,criterion=None):
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        sub_arrays = JunctionCountMatrixQuery.split_array_to_chunks(self.results[0],self.cores)
        sub_conds = JunctionCountMatrixQuery.split_df_to_chunks(self.cond_subset_df,self.cores)
        hlas = self.results[1]
        pool = mp.Pool(processes=self.cores)

        r = [pool.apply_async(func=JunctionCountMatrixQuery.show_neoantigen_burden_single_run,args=(sub_array,sub_cond,hlas,stage,verbosity,contain_uid,criterion,)) for sub_array,sub_cond in zip(sub_arrays,sub_conds)]
        pool.close()
        pool.join()
        results = []
        for collect in r:
            result = collect.get()
            results.append(result)
        
        burden_valid = pd.concat(results,axis=0)
        burden_valid.index = self.subset.index
        burden_valid.columns = self.subset.columns

        burden_invalid = self.junction_count_matrix.loc[self.invalid,:].copy()
        for col in burden_invalid.columns:
            burden_invalid.loc[:,col] = np.full(burden_invalid.shape[0],0)
        
        burden = pd.concat([burden_valid,burden_invalid],axis=0)

        # statistics
        burden['mean'] = burden.mean(axis=1)
        burden.loc['burden'] = burden.sum(axis=0)
        burden.sort_values(by='mean',axis=0,ascending=False,inplace=True)
        burden.sort_values(by='burden',axis=1,ascending=True,inplace=True)
        # reorder a bit
        c_list = burden.columns.tolist()
        c_list.remove('mean')
        c_list.append('mean')
        i_list = burden.index.tolist()
        i_list.remove('burden')
        i_list.append('burden')
        burden = burden.loc[i_list,c_list]
        burden.to_csv(os.path.join(outdir,name),sep='\t')    

    def show_neoantigen_frequency(self,outdir,name,stage,verbosity,contain_uid,plot,plot_name=None,yscale='linear',criterion=None):
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        sub_arrays = JunctionCountMatrixQuery.split_array_to_chunks(self.results[0],self.cores)
        sub_conds = JunctionCountMatrixQuery.split_df_to_chunks(self.cond_subset_df,self.cores)
        hlas = self.results[1]
        column_names = self.subset.columns
        pool = mp.Pool(processes=self.cores)

        r = [pool.apply_async(func=JunctionCountMatrixQuery.show_neoantigen_frequency_single_run,args=(sub_array,sub_cond,hlas,column_names,stage,verbosity,contain_uid,criterion,)) for sub_array,sub_cond in zip(sub_arrays,sub_conds)]
        pool.close()
        pool.join()
        results = []
        for collect in r:
            result = collect.get()
            results.append(result)  

        # let's define another function for merging dictionary, but not overwrite the values
        def merge_key_and_values(dic1,dic2):
            return {x:list(set(dic1.get(x,[]) + dic2.get(x,[]))) for x in set(dic1).union(set(dic2))}         

        dic = reduce(merge_key_and_values,results) 

        # statistics
        df = pd.Series(data=dic).to_frame(name='samples')
        df['n_sample'] = df.apply(lambda x:len(set(x[0])),axis=1).values
        df.sort_values(by='n_sample',ascending=False,inplace=True)
        df.to_csv(os.path.join(outdir,name),sep='\t')
        # plot
        if plot:
            fig,ax = plt.subplots()
            ax.bar(x=np.arange(df.shape[0]),height=df['n_sample'].values,edgecolor='k')
            ax.set_xlabel('Neoantigen rank by its occurence (descending order)')
            ax.set_ylabel('Occurence (n_sample)')
            ax.set_title('Neoantigen Occurence')
            plt.savefig(os.path.join(outdir,'x_neoantigen_{}'.format(plot_name)),bbox_inches='tight')
            plt.close()
            fig,ax = plt.subplots()
            try:
                sns.histplot(df['n_sample'].values,binwidth=1,kde=True,ax=ax)
                ax.set_yscale(yscale)
                plt.savefig(os.path.join(outdir,'x_occurence_{}'.format(plot_name)),bbox_inches='tight')
                plt.close()
            except ValueError:
                print('All neoantigens are present in same amount of patients, seaborn histplot can not properly figure out binedge, no occurence plot generated')





    def show_neoantigen_as_fasta(self,outdir,name,stage,verbosity,contain_uid,sample=None,criterion=None):
        '''
        write the neoantigen as a fasta file for MS validation or other purpose

        :param outdir: string ,the output directory for the generated fasta file
        :param name: string, the name of the output fasta file
        :param stage: int, either 0,1,2,3, the stage of neoantigen you want to report
        :param verbosity: int, either 1,2,3, the verbosity of candidates
        :param contain_uid: boolean, whether you want to contain the junction UID along with the neoantigen
        :param sample: string, the name of the sample in which you want to extract the neoantigen
        :param criterion: nested list, the criterion for filtering out the candidate

        Examples::

            jcmq = snaf.JunctionCountMatrixQuery.deserialize('result/after_prediction.p')
            sample = 'SRR5933735.Aligned.sortedByCoord.out'
            jcmq.show_neoantigen_as_fasta(outdir='./fasta',name='neoantigen_{}.fasta'.format(sample),stage=2,verbosity=1,contain_uid=True,sample=sample)

        '''
        if not os.path.exists(outdir):
            os.mkdir(outdir)
            
        sub_arrays = JunctionCountMatrixQuery.split_array_to_chunks(self.results[0],self.cores)
        sub_conds = JunctionCountMatrixQuery.split_df_to_chunks(self.cond_subset_df,self.cores)
        hlas = self.results[1]
        column_names = self.subset.columns
        pool = mp.Pool(processes=self.cores)

        r = [pool.apply_async(func=JunctionCountMatrixQuery.show_neoantigen_frequency_single_run,args=(sub_array,sub_cond,hlas,column_names,stage,verbosity,contain_uid,criterion,)) for sub_array,sub_cond in zip(sub_arrays,sub_conds)]
        pool.close()
        pool.join()
        results = []
        for collect in r:
            result = collect.get()
            results.append(result)  

        # let's define another function for merging dictionary, but not overwrite the values
        def merge_key_and_values(dic1,dic2):
            return {x:list(set(dic1.get(x,[]) + dic2.get(x,[]))) for x in set(dic1).union(set(dic2))}         

        dic = reduce(merge_key_and_values,results) 
        # statistics
        df = pd.Series(data=dic).to_frame(name='samples')
        df['n_sample'] = df.apply(lambda x:len(set(x[0])),axis=1).values
        df.sort_values(by='n_sample',ascending=False,inplace=True)
        # subset
        if sample is not None:
            df = df.loc[df.apply(lambda x:sample in x['samples'], axis=1),:] 
            name = 'neoantigen_' + sample + '.fasta'
        # write to fasta
        with open(os.path.join(outdir,name),'w') as f:
            for row in df.itertuples():
                if not contain_uid:   # as the row.Index is a tuple
                    tmp = [str(item) for item in row.Index]
                    f.write('>{}|{}\n'.format('|'.join(tmp),row.n_sample)) 
                    f.write('{}\n'.format(row.Index[0]))   
                else:   # as the row.Index is a string delimited by comma
                    f.write('>{}|{}\n'.format(row.Index.replace(',','|'),row.n_sample)) 
                    f.write('{}\n'.format(row.Index.split(',')[0]))
                


    # let's define an atomic function inside here
    @staticmethod
    def show_neoantigen_burden_single_run(sub_array,sub_cond,hlas,stage,verbosity,contain_uid,criterion=None):
        nj_burden_2d = []
        for nj,index in zip(sub_array,range(sub_cond.shape[0])):
            nj_burden = []
            cond_row = sub_cond.iloc[index].values
            for hla,cond in zip(hlas,cond_row):
                if nj is None or not cond :
                    nj_burden.append(0)
                else:
                    nj_copy = deepcopy(nj)
                    nj_copy.enhanced_peptides = nj_copy.enhanced_peptides.filter_based_on_hla(selected_hla=hla)
                    nj_copy.derive_candidates(stage=stage,verbosity=verbosity,contain_uid=contain_uid,criterion=criterion)
                    nj_burden.append(len(nj_copy.candidates))
            nj_burden_2d.append(nj_burden)  
        df = pd.DataFrame(data=nj_burden_2d)
        return df

    # let's define the single function to run in each subprocess
    @staticmethod
    def show_neoantigen_frequency_single_run(sub_array,sub_cond,hlas,column_names,stage,verbosity,contain_uid,criterion=None):
        dic = {}
        for nj,index in zip(sub_array,range(sub_cond.shape[0])):
            cond_row = sub_cond.iloc[index].values
            for hla,column_name,cond in zip(hlas,column_names,cond_row):
                if nj is not None and cond:
                    nj_copy = deepcopy(nj)
                    nj_copy.enhanced_peptides = nj_copy.enhanced_peptides.filter_based_on_hla(selected_hla=hla)
                    nj_copy.derive_candidates(stage=stage,verbosity=verbosity,contain_uid=contain_uid,criterion=criterion)
                    for cand in nj_copy.candidates:
                        try:
                            dic[cand].append(column_name)
                        except KeyError:
                            dic[cand] = []
                            dic[cand].append(column_name)       
        return dic   

class EnhancedPeptides():
    def __init__(self,peptides,hlas,kind):
        # kind means how to instantiate the EnhancedPeptides object
        if kind == 2:  # hla-reduced enhanced peptide, this mode is used when filtsre_based_on_hla
            self.mers = peptides
            self.info = hlas
        else:
            self.mers = []
            self.info = []
            for k,v in peptides.items():
                self.mers.append(k)
                phlas = {}  # {pep1:{origin:(extra,n_from_first,phase,evidences),hla1:{},hla2:{}}}
                for pep in v:
                    if kind == 0:  # each pep is (pep,extra,n_from_first,phase,evidences), this mode is used before binding_prediction
                        pairs = {}
                        pairs['origin'] = (pep[1],pep[2],pep[3],pep[4])
                        for hla in hlas:
                            pairs[hla] = {}
                        phlas[pep[0]] = pairs
                    elif kind == 1:   # each pep is (pep,extra,n_from_first,hla,phase,evidences), this mode is used when filter_based_on_criterion
                        try:
                            phlas[pep[0]]['origin'] = (pep[1],pep[2],pep[4],pep[5])
                        except KeyError:
                            phlas[pep[0]] = {}
                            phlas[pep[0]]['origin'] = (pep[1],pep[2],pep[4],pep[5])
                        phlas[pep[0]][pep[3]] = {}
                self.info.append(phlas)

    def __str__(self):
        return str(self.info)

    def __getitem__(self,key):
        index = self.mers.index(key)
        return self.info[index]

    def is_empty(self):
        total = 0
        for dic in self.info:
            total += len(dic)
        if total == 0:
            return True
        else:
            return False

    def simplify_to_list(self,verbosity):   # make sure to first filter and re-instantiate
        if verbosity==1:  # only peptide, result is also tuple
            result_list = []
            for mer in self.mers:
                result_list.extend(list(self[mer].keys()))   # utlize __getitem__ method, so self[9] instead of self[0]
            result_list = [(item,) for item in result_list]
        elif verbosity==2:  # peptide and hla
            result_list = []
            for mer in self.mers:
                for pep,hla_complex in self[mer].items():
                    for hla in hla_complex.keys():
                        if hla != 'origin':
                            result_list.append((pep,hla))
        elif verbosity==3:
            result_list = []
            for mer in self.mers:
                for pep,hla_complex in self[mer].items():
                    extra,n_from_first,phase,evidences = hla_complex['origin']
                    for hla in hla_complex.keys():
                        if hla != 'origin':
                            result_list.append((pep,hla,extra,n_from_first))
        return result_list

    def register_attr(self,df,attr_name):
        '''
        df should follow:
          peptide     mer     hla       score    identity
         AAAAAAAAA     9   HLA-A*01:01   0.3       SB        
        '''
        for mer,sub_df in df.groupby(by='mer'):
            index = self.mers.index(int(mer))
            for peptide,sub_sub_df in sub_df.groupby(by='peptide'):
                for row in sub_sub_df.itertuples(index=False):
                    self.info[index][peptide][row.hla][attr_name] = (float(row.score),str(row.identity))

    def filter_based_on_hla(self,selected_hla,reinstantiate=True):
        new_mers = []
        new_info = []
        selected_hla = set(selected_hla)
        for i,k in enumerate(self.mers):
            new_info_container = {}
            for pep,hla_complex in self.info[i].items():
                reduced_hla_complex = {}
                for hla,attrs in hla_complex.items():
                    if hla == 'origin':
                        reduced_hla_complex['origin'] = attrs
                    else:
                        if hla in selected_hla:
                            reduced_hla_complex[hla] = attrs
                new_info_container[pep] = reduced_hla_complex
            new_mers.append(k)
            new_info.append(new_info_container)
        if reinstantiate:
            ep = EnhancedPeptides(peptides=new_mers,hlas=new_info,kind=2)
        return ep
                

        


    def filter_based_on_criterion(self,criteria,reinstantiate=True):
        # criterion: [(net),], ['netMHCpan_el',1,==,SB]
        peptides = {k:[] for k in self.mers}
        for i,k in enumerate(self.mers):
            for pep,hla_complex in self.info[i].items():
                extra,n_from_first,phase,evidences = hla_complex['origin']
                for hla,attrs in hla_complex.items():
                    if hla == 'origin':
                        continue
                    boolean_list = []
                    for criterion in criteria:
                        if criterion[1] == 1:   # matched is not a number
                            eval_string = 'attrs[\'{}\'][{}] {} \'{}\''.format(criterion[0],criterion[1],criterion[2],criterion[3])
                        elif criterion[1] == 0:   # matched is a number
                            eval_string = 'attrs[\'{}\'][{}] {} {}'.format(criterion[0],criterion[1],criterion[2],criterion[3])
                        try:
                            boolean = eval(eval_string)                            
                        except KeyError:
                            boolean = False
                        boolean_list.append(boolean)
                    boolean_final = all(boolean_list)
                    if boolean_final:
                        peptides[k].append((pep,extra,n_from_first,hla,phase,evidences))   # if not reinstantiate, this format is called reduced form
        if reinstantiate:
            return EnhancedPeptides(peptides,None,1)
        else:
            return peptides


class NeoJunction():
    def __init__(self,uid,count,check_gtex):
        if check_gtex:
            NeoJunction.is_neojunction(uid,count)
        self.uid = uid
        self.count = count

    def __str__(self):
        try:
            ts = self.tumor_specificity
        except AttributeError:
            ts = None
        try:
            et = self.event_type
        except AttributeError:
            et = None
        try:
            j = self.junction[:5] + '...' + self.junction[-5:]
        except AttributeError:
            j = None
        try:
            peptides = {k:len(v) for k,v in self.peptides.items()}
        except AttributeError:
            peptides = None
        try:
            ep = {i:len(self.enhanced_peptides[i]) for i in self.enhanced_peptides.mers}
        except AttributeError:
            ep = None
        try:
            cand_len = len(self.candidates)
            cand_example = self.candidates[0]
        except AttributeError:
            cand_len,cand_example = None,None
        return 'uid: {}\n'\
               'count: {}\n'\
               'tumor specificity: {}\n'\
               'event type: {}\n'\
               'junction: {}\n'\
               'peptides: {}\n'\
               'Enhanced peptides: {}\n'\
               'Candidates: len is {}, first is {}\n'.format(self.uid,self.count,ts,et,j,peptides,ep,cand_len,cand_example)



    @staticmethod
    def is_neojunction(uid,count):
        identity,detail = crude_tumor_specificity(uid,count)
        if not identity:
            raise Exception('This is not a NeoJunction, instantiation fails')

    def infer_tumor_specificity(self,method):
        tumor_specificity = accurate_tumor_specificity(self.uid,method)
        self.tumor_specificity = tumor_specificity
        return tumor_specificity

    def gtex_viewer(self,kind,outdir):
        if kind == 1:   # line plot, all tissues, count (norm or not)
            gtex_visual_count(self.uid,norm=True,out_folder=outdir)
        elif kind == 2: # hist plot, combined, norm count
            gtex_visual_norm_count_combined(self.uid,out_folder=outdir)
        elif kind == 3:  # hist plot, per tissue, number of samples that are non-zero in each tissue
            gtex_visual_per_tissue_count(self.uid,out_folder=outdir)

    
    def detect_type(self):
        '''
        Ordinary: ENSG00000107902:E10.1-E12.1
        Alt3: ENSG00000110057:E5.1-E6.2_67996641
        Alt5: ENSG00000100321:E7.1_39364266-E8.1
        Intron Retention: ENSG00000115524:I4.1-E5.1
        Novel Exon: ENSG00000008441:I40.1_13076665-E41.1
        Trans-splicing: ENSG00000196565:E14.2-ENSG00000213934:E3.1
        UTR Event: ENSG00000164068:U0.1_49689185-E2.1
        '''
        valid_pattern = re.compile(r'^ENSG\d+:.+?-.+')
        if re.search(valid_pattern,self.uid):   # at least valid one
            if len(re.findall('ENSG',self.uid)) == 2:
                event_type = 'trans_splicing'
            elif 'U' in self.uid:
                event_type = 'utr_event'
            elif '_' in self.uid:
                subexon12 = self.uid.split(':')[1]
                subexon1, subexon2 = subexon12.split('-')
                if 'I' in subexon12:
                    event_type = 'novel_exon'
                elif '_' in subexon1 and '_' in subexon2:
                    event_type = 'alt5_alt3'
                elif '_' in subexon1 and '_' not in subexon2:
                    event_type = 'alt5'
                elif '_' in subexon2 and '_' not in subexon1:
                    event_type = 'alt3'
                else:
                    event_type = 'invalid'
            elif 'I' in self.uid:
                event_type = 'intron_retention'
            elif re.search(r'^ENSG\d+:E\d+\.\d+-E\d+\.\d+$',self.uid):
                e = self.uid.split(':')[1]
                e1 = e.split('-')[0]
                e2 = e.split('-')[1]
                e1_int = int(e1.split('.')[0][1:])
                e1_frac = int(e1.split('.')[1])
                e2_int = int(e2.split('.')[0][1:])
                e2_frac = int(e2.split('.')[1])                
                if e1 == e2:   # E5.1-E5.1
                    event_type = 'invalid'
                else:
                    if e1_int > e2_int or (e1_int==e2_int and e1_frac>e2_frac):   # E13.1-E12.4   or E12.10-E12.9
                        event_type = 'invalid'
                    else:
                        event_type = 'ordinary'
            else:
                event_type = 'invalid'
        else:
            event_type = 'invalid'
        self.event_type = event_type
        return event_type

    def retrieve_junction_seq(self):
        if self.event_type != 'invalid':
            ensid = self.uid.split(':')[0]
            subexon1,subexon2 = ':'.join(self.uid.split(':')[1:]).split('-')
            code = is_consecutive(subexon1,subexon2)
            seq1 = subexon_tran(subexon1,ensid,'site1',code)
            seq2 = subexon_tran(subexon2,ensid,'site2',code)
            junction = ','.join([seq1,seq2])
            self.junction = junction
        else:
            self.junction = '$' * 10   # indicating invalid uid


    def in_silico_translation(self,ks=[9,10],strict=False):
        peptides = {k:[] for k in ks}
        coord = uid_to_coord(self.uid)  # chr1:555-666(+)
        if '$' not in self.junction and '*' not in self.junction and '#' not in self.junction and 'unknown' not in coord:
            first,second = self.junction.split(',')
            # annotate peptide based on starting codon evidence
            ensg = self.uid.split(':')[0]
            strand = coord.split('(')[1].rstrip(')')
            if strand == '+':
                coord_first_exon_last_base = coord.split(':')[1].split('-')[0]
            elif strand == '-':
                coord_first_exon_last_base = coord.split(':')[1].split('(')[0].split('-')[1]
            possible_start_codon_coord = dict_start_codon.get(ensg,[])   # [333,444] or []
            support_phases_dict = {}   # {0:[(333,enst,+)]} or {}
            for pssc in possible_start_codon_coord:
                supports = get_support_phase(ensg,coord_first_exon_last_base,pssc,strand,len(first))
                for (phase_, pssc_, enst_, strand_) in supports:
                    support_phases_dict.setdefault(phase_,[]).append((pssc_,enst_,strand_))
            for phase in [0,1,2]:  # tranlation starts with index "phase"
                evidences = tuple(support_phases_dict.get(phase,[]))
                if strict and len(evidences)==0:  # in strict mode, without start_codon evidence, let's skip this phase
                    continue
                de_facto_first = first[phase:]
                pep_dict = get_peptides(de_facto_first,second,ks,phase,evidences)   # (pep,extra,n_from_first,phase,evidences)
                for k in ks:  
                    peptides[k].extend(pep_dict[k])
                    peptides[k] = list(set(peptides[k]))
        self.peptides = peptides
        return peptides

    def binding_prediction(self,hlas,binding_method=None):
        hlas = list(set(hlas))
        ep = EnhancedPeptides(self.peptides,hlas,0)
        if ep.is_empty():
            raise Exception('Already no candidates after in-silico translation')
        for k,v in self.peptides.items():
            try:
                v = list(zip(*v))[0]
            except IndexError:
                continue
            if binding_method == 'netMHCpan':
                df = run_netMHCpan(software_path,v,hla_formatting(hlas,'netMHCpan_output','netMHCpan_input'),k)
            elif binding_method == 'MHCflurry':
                df = run_MHCflurry(v,hla_formatting(hlas,'netMHCpan_output','deepimmuno'))
                df['hla'] = hla_formatting(df['hla'].tolist(),'deepimmuno','netMHCpan_output')
            ep.register_attr(df,'netMHCpan_el')
        self.enhanced_peptides = ep


    def immunogenicity_prediction(self):
        reduced = self.enhanced_peptides.filter_based_on_criterion([('netMHCpan_el',0,'<=',2),],False)
        ep = self.enhanced_peptides.filter_based_on_criterion([('netMHCpan_el',0,'<=',2),],True)
        if ep.is_empty():
            raise Exception('Already no candidates after binding prediction')
        for k,v in reduced.items():
            try:
                v_pep,v_hla = list(zip(*v))[0],list(zip(*v))[3]
            except IndexError:  # one of the mer, no candicate, but the other mer has, so pass the empty check
                continue
            data = np.column_stack((v_pep,v_hla))
            df_input = pd.DataFrame(data=data)
            df_input[1] = hla_formatting(df_input[1].tolist(),'netMHCpan_output','deepimmuno')
            df_output = run_deepimmuno(df_input)
            df = pd.DataFrame({'peptide':df_output['peptide'].values,'mer':[k]*df_output.shape[0],
                               'hla':hla_formatting(df_output['HLA'].values.tolist(),'deepimmuno','netMHCpan_output'),
                               'score':df_output['immunogenicity'].values,
                               'identity':[True if item > 0.5 else False for item in df_output['immunogenicity'].values]})
            self.enhanced_peptides.register_attr(df,attr_name='deepimmuno_immunogenicity')


    def derive_candidates(self,stage,verbosity,contain_uid=True,criterion=None):
        if stage == 1: # all translated peptides
            self.candidates = self.enhanced_peptides.simplify_to_list(verbosity=verbosity)
        elif stage == 2: # all bound peptides
            self.candidates = self.enhanced_peptides.filter_based_on_criterion([('netMHCpan_el',0,'<=',2),]).simplify_to_list(verbosity=verbosity)
        elif stage == 3: # immunogenic peptides
            self.candidates = self.enhanced_peptides.filter_based_on_criterion([('netMHCpan_el',0,'<=',2),('deepimmuno_immunogenicity',1,'==','True'),]).simplify_to_list(verbosity=verbosity)
        elif stage == 'custom':
            self.candidates = self.enhanced_peptides.filter_based_on_criterion(criterion).simplify_to_list(verbosity=verbosity)
        if contain_uid:
            new = []
            for item in self.candidates:
                tmp = [str(i) for i in item]
                tmp.append(self.uid)
                new.append(','.join(tmp))
            self.candidates = new





    def visualize(self,outdir,name,criterion=[('netMHCpan_el',0,'<=',2),('deepimmuno_immunogenicity',1,'==','True'),]):
        reduced = self.enhanced_peptides.filter_based_on_criterion(criterion,False) # now each element should add phase and evidences as well
        '''
        {9: [('LPSPPAQEL', 2, 0, 'HLA-B*08:01'), ('LPSPPAQEL', 2, 0, 'HLA-B*08:02'), 
             ('SLYLLLQHR', 1, 2, 'HLA-A*68:01')], 
        10: [('TSLYLLLQHR', 1, 3, 'HLA-A*68:01'), 
             ('TLPSPPAQEL', 2, 1, 'HLA-A*02:01'), ('TLPSPPAQEL', 2, 1, 'HLA-A*24:02'), ('TLPSPPAQEL', 2, 1, 'HLA-B*08:01'), 
             ('TLPSPPAQEL', 2, 1, 'HLA-B*08:02')]}
        '''
        n_axes = 0
        to_draw = []
        for k,v in reduced.items():
            to_draw.extend(v)
            n_axes += len(v)
        ncols = 4
        nrows = n_axes // ncols + 1 + 1
        height_ratios = [0.2 if i == 0 else (1-0.2)/(nrows-1) for i in range(nrows)]
        fig = plt.figure()
        gs = mpl.gridspec.GridSpec(nrows=nrows,ncols=ncols,height_ratios=height_ratios,wspace=0.5,hspace=0.5)
        ax_genome  = fig.add_subplot(gs[0,:])
        axes_list = []
        for i in range(1,gs.nrows):
            for j in range(0,gs.ncols):
                ax = fig.add_subplot(gs[i,j])
                axes_list.append(ax)
        ax_genome = draw_genome(ax_genome,self.uid,dict_exonCoords)
        for i,ax in enumerate(axes_list):
            if i < n_axes:
                # info
                aa, extra, n_from_first, hla, phase, evidences = to_draw[i]
                first,second = self.junction.split(',')
                dna_first = extra + n_from_first * 3
                dna_second = -extra + (len(aa)-n_from_first) * 3
                binding_score = self.enhanced_peptides[len(aa)][aa][hla]['netMHCpan_el'][0]
                immunogenicity_score = self.enhanced_peptides[len(aa)][aa][hla]['deepimmuno_immunogenicity'][0]   
                # draw     
                ax = show_candicates(ax,aa,extra,n_from_first,hla,phase,evidences,first,second,dna_first,dna_second,binding_score,immunogenicity_score) 
            else:
                ax.axis('off') 
        fig.subplots_adjust(top=0.9)             
        fig.suptitle('{} Count:{}'.format(self.uid,self.count))
        plt.savefig(os.path.join(outdir,name),bbox_inches='tight')
        plt.close()

       

# processing functions
def hla_formatting(pre,pre_type,post_type):
    if pre_type == 'netMHCpan_output' and post_type == 'netMHCpan_input':  # HLA-A*01:01 to HLA-A01:01
        post = [hla.replace('*','') for hla in pre]
    elif pre_type == 'netMHCpan_input' and post_type == 'netMHCpan_output':  # HLA-A01:01 to HLA-A*01:01
        post = [hla[:5] + '*' + hla[5:] for hla in pre]
    elif pre_type == 'netMHCpan_output' and post_type == 'deepimmuno':  # HLA-A*01:01 to HLA-A*0101
        post = [hla.replace(':','') for hla in pre]
    elif pre_type == 'deepimmuno' and post_type == 'netMHCpan_output': # HLA-A*0101 to HLA-A*01:01
        post = [hla[:8] + ':' + hla[-2:] for hla in pre]
    return post





def get_peptides(de_facto_first,second,ks,phase,evidences):
    peptides = {k:[] for k in ks}
    extra = len(de_facto_first) % 3  # how many base left in first assuming no stop condon in front of it.
    num = len(de_facto_first) // 3   # how many set of codons in the first.
    aa_first = str(Seq(de_facto_first).translate(to_stop=True))  # if not the multiple of 3, only translate the part that are divisible by 3, it will be a warning for BioPython
    if len(aa_first) == num:  # successfully read through
        if extra == 0:
            continue_second = second
        elif extra == 1:
            continue_second = de_facto_first[-1] + second
        elif extra == 2:
            continue_second = de_facto_first[-2:] + second
        aa_second = str(Seq(continue_second).translate(to_stop=True))
        if len(aa_second) > 0:  # at least, not ''
            for k in ks:
                second_most = min(k,len(aa_second)) # the max allowed number of residue for aa_second
                first_most = len(aa_first)  # the max allowed number of residue for aa_first
                for n_from_second in range(second_most,0,-1):
                    n_from_first = k - n_from_second
                    if n_from_first == 0 and extra == 0:
                        '''
                        cttca cct cac ttt acc ttc tcc tcc agc aca gga act agg aac tac gga gag aga agc caa 
                           S   P   H   F   T   F   S   S   S   T   G   T   R   N   Y   G   E   R   S   Q

                        the comma is between "ttt" and "acc"
                        so, when extra == 0, means whole second will be after the junction, in this case, we have to require 
                        the peptide at least include a amino acid from first part, otherwise, TFSSSTGTR won't be a splice-derived
                        peptide.   

                        this example can be reproduced using: nj = NeoJunction('ENSG00000223572:E15.1-E15.2')                     
                        '''
                        continue
                    if n_from_first <= first_most:
                        if n_from_first > 0:
                            pep = aa_first[-n_from_first:] + aa_second[:n_from_second]
                        elif n_from_first == 0:
                            pep = aa_second[:n_from_second]
                        peptides[k].append((pep,extra,n_from_first,phase,evidences))
    return peptides
                        


def query_from_dict_fa(abs_start,abs_end,EnsID,strand):
    '''
    abs_start and abs_end always means the xth base in forward strand

    the returned exon_seq, however, means the 5'-3' seq depending on the strand information.
    '''
    if strand == '+':        
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq = dict_fa[EnsID][3]
        start_index = int(abs_start) - start + 2000
        end_index = int(abs_end) - start + 1 + 2000
        exon_seq = seq[start_index:end_index]
    
    elif strand == '-':
        start = int(dict_fa[EnsID][1])
        end = int(dict_fa[EnsID][2])
        seq_reverse = dict_fa[EnsID][3]
        seq_forward = str(Seq(seq_reverse).reverse_complement())  # Hs_gene.fa restore the reverse strand info
        start_index = int(abs_start) - start + 2000
        end_index = int(abs_end) - start + 1 + 2000 # endpoint in python is non-inclusive
        exon_seq_1 = seq_forward[start_index:end_index]
        s = Seq(exon_seq_1)
        exon_seq = str(s.reverse_complement())
    return exon_seq

            


def utrAttrs(EnsID):  # try to get U0.1's attribute, but dict_exonCoords doesn't have, so we just wanna get the first entry for its EnsGID
    exonDict = dict_exonCoords[EnsID] 
    attrs = next(iter(exonDict.values()))
    chr_,strand = attrs[0],attrs[1]
    return chr_,strand

def utrJunction(site,EnsGID,strand,chr_,flag,seq_len=100):  # U0.1_438493849, here 438493849 means the site (suffix)
    if flag == 'site1' and strand == '+':  # U0.1_438493849 - E2.2
        otherSite = int(site) - seq_len + 1   # extract UTR with length = 100
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(otherSite),int(site))
    elif flag == 'site1' and strand == '-':    
        otherSite = int(site) + seq_len - 1 
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(site),int(otherSite))
        exon_seq = str(Seq(exon_seq).reverse_complement())
    elif flag == 'site2' and strand == '+':  # E5.3 - U5.4_48374838
        otherSite = int(site) + seq_len -1
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(site),int(otherSite))
    elif flag == 'site2' and strand == '-':
        otherSite = int(site) - seq_len + 1
        exon_seq = retrieveSeqFromUCSCapi(chr_,int(otherSite),int(site))
        exon_seq = str(Seq(exon_seq).reverse_complement())
    return exon_seq

def retrieveSeqFromUCSCapi(chr_,start,end):
    url = 'http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment={0}:{1},{2}'.format(chr_,start,end)
    response = requests.get(url)
    status_code = response.status_code
    assert status_code == 200
    try:
        my_dict = xmltodict.parse(response.content)
    except:
        exon_seq = '#' * 10  # indicating the UCSC doesn't work
        return exon_seq
    exon_seq = my_dict['DASDNA']['SEQUENCE']['DNA']['#text'].replace('\n','').upper()
    return exon_seq

def is_consecutive(subexon1,subexon2):
    # I only care of E3.3-E3.4, other form will automatically being labelled as not consecutive
    # the return code: 0 -> not consecutive, 1 -> consecutive
    code = 0
    pattern1 = re.compile(r'^E(\d{1,3})\.(\d{1,3})$')
    match1 = re.search(pattern1,subexon1)
    if match1:
        captured_group1 = match1.group(1)
        captured_group2 = match1.group(2)
        pattern2 = re.compile(rf'^E{captured_group1}\.{int(captured_group2)+1}$')
        match2 = re.search(pattern2,subexon2)
        if match2:
            code = 1
    return code

def subexon_tran(subexon,EnsID,flag,code):  # flag either site1 or site2
    '''
    1. subexon can take multiple forms depending on the event type
    E1.2 or I3.4
    E6.2_67996641 or I40.1_13076665, also depending on whether they are subexon1 or subexon2
    ENSG00000213934:E3.1 or ENSG00000213934:E2.1_473843893894
    U0.1_49689185
    2. everything with trailing suffix will depend on the subexon1 or subexon2, but sometimes, it is fixed (trans-splicing can only be in subexon2)
    3. to be clear, the exon_seq returned is always 5'-3' sequence, not forward anymore.

    if code is 1, meaning its consecutive, E3.3-E3.4, if the flag is site2, then offset the start position of it by 1, the start position is the conceptual
    start position, if it is in negative string, it will be operated by the stop position.
    '''
    # let's remove those in the first place
    # if EnsID not in set(dict_exonCoords.keys()):
    #     exon_seq = '*' * 10  # indicator for error on MultiPath-PSI itself
    #     return exon_seq
        
    try:   # E1.2 or I3.4
        attrs = dict_exonCoords[EnsID][subexon]  # [chr,strand,start,end,suffer] 
        if code == 1 and flag == 'site2':
            if attrs[1] == '+':
                exon_seq = query_from_dict_fa(int(attrs[2])+1,attrs[3],EnsID,attrs[1])   
            else:
                exon_seq = query_from_dict_fa(attrs[2],int(attrs[3])-1,EnsID,attrs[1])
        else:
            exon_seq = query_from_dict_fa(attrs[2],attrs[3],EnsID,attrs[1])
    except KeyError:
        if ':' in subexon: # ENSG00000213934:E3.1
            fusionGeneEnsID = subexon.split(':')[0] 
            fusionGeneExon = subexon.split(':')[1]        
            if  '_' in fusionGeneExon:   # ENSG:E2.1_473843893894
                suffix = fusionGeneExon.split('_')[1]
                subexon = fusionGeneExon.split('_')[0]
                attrs = dict_exonCoords[fusionGeneEnsID][subexon]
                if attrs[1] == '+':  
                    exon_seq = query_from_dict_fa(suffix,attrs[3],fusionGeneEnsID,attrs[1]) 
                else:  
                    exon_seq = query_from_dict_fa(attrs[2],suffix,fusionGeneEnsID,attrs[1])
            else:  # ENSG:E2.1
                try:
                    attrs = dict_exonCoords[fusionGeneEnsID][fusionGeneExon]
                except KeyError:
                    exon_seq = '*' * 10  # indicator for error on MultiPath-PSI itself
                else:
                    exon_seq = query_from_dict_fa(attrs[2],attrs[3],fusionGeneEnsID,attrs[1]) 

        else:  # could be trailing or utr, or non-existing ordinary subexon
            try:
                suffix = subexon.split('_')[1]
            except IndexError: # the logic is there's a subexon E45.3, it is no trailing, but just not in the exonCoords.
                exon_seq = '*' * 10  # indicator for error on MultiPath-PSI itself
            else:
                subexon = subexon.split('_')[0]
                try:
                    attrs = dict_exonCoords[EnsID][subexon]
                except KeyError:  # must be UTR
                    chrUTR,strandUTR = utrAttrs(EnsID) # this is get from a random subexon under that EnsID
                    exon_seq = utrJunction(suffix,EnsID,strandUTR,chrUTR,flag)  
                else:   # must be trailing
                    if flag == 'site2':
                        if attrs[1] == '+':  
                            exon_seq = query_from_dict_fa(suffix,attrs[3],EnsID,attrs[1]) 
                        else:  
                            exon_seq = query_from_dict_fa(attrs[2],suffix,EnsID,attrs[1])
                    elif flag == 'site1':  # not affected by overhang since it is site1
                        if attrs[1] == '+': 
                            exon_seq = query_from_dict_fa(attrs[2],suffix,EnsID,attrs[1])
                        else:
                            exon_seq = query_from_dict_fa(suffix,attrs[3],EnsID,attrs[1])
    return exon_seq


def enhance_frequency_table(df,remove_quote=True,save=True,outdir='',name=None):
    '''
    This is a wrapper function to add (1) gene symbol, (2) chromosome coordinate, (3) tumor specificity to the frequency table at once

    :param df: the dataframe for the frequency table (freq_stage{n}_verbosity1_uid.txt)
    :param remove_quote: boolean, depending on how your df is loaded, if directly read from the disk, certain column will contain quotation, whether to remove it or not
    :param save: boolean, whether to save the result or not, default is True
    :param outdir: string, the output folder.
    :param name: string, the output file name

    Example::

        enhance_frequency_table(df=df,True,True,'result','freq_stage3_verbosity1_uid_gene_symbol_coord_mean_mle.txt')

    '''
    print('adding gene symbol')
    df = add_gene_symbol_frequency_table(df=df,remove_quote=remove_quote)
    print('adding chromosome coordinates')
    df = add_coord_frequency_table(df,remove_quote=False)
    print('adding tumor specificity mean raw count')
    df = add_tumor_specificity_frequency_table(df,'mean',False)
    print('adding tumor specificity MLE score')
    df = add_tumor_specificity_frequency_table(df,'mle',False)
    if save:
        df.to_csv(os.path.join(outdir,name),sep='\t')
    return df

def add_coord_frequency_table(df,remove_quote=True):
    '''
    Convert the uid to chromsome coordinates

    :param df: the input df, index must be uid
    :param remove_quote: bool, depending on how your df is loaded, if directly read from the disk, certain column will contain quotation, whether to remove it or not
    :return df: the output df, with added coordinate column

    Exmaple::

        add_coord_frequency_table(df=df,remove_quote=True)
    '''
    # the index has to be comma separated string
    from ast import literal_eval
    if remove_quote:
        df['samples'] = [literal_eval(item) for item in df['samples']]   
    uid_list = [item.split(',')[1] for item in df.index] 
    coord_list = [uid_to_coord(uid) for uid in uid_list]
    df['coord'] = coord_list
    return df


def uid_to_coord(uid):
    tmp_list = uid.split(':')
    if len(tmp_list) == 2:
        ensg,exons = tmp_list
    elif len(tmp_list) == 3:
        ensg = tmp_list[0]
        exons = ':'.join(tmp_list[1:])
    first,second = exons.split('-')
    # figure out start_coord
    if '_' in first:
        actual_exon,trailing = first.split('_')
        try:
            attrs = dict_exonCoords[ensg][actual_exon]
        except KeyError:
            if 'U' in actual_exon:
                proxy_exon = list(dict_exonCoords[ensg].keys())[0]
                attrs = dict_exonCoords[ensg][proxy_exon]
                chrom = attrs[0]
                strand = attrs[1]
                start_coord = trailing
            else:   # probably a rare error
                chrom = 'unknown'
                strand = 'unknown'
                start_coord = 'unknown'
        else:
            chrom = attrs[0]
            strand = attrs[1]
            start_coord = trailing
    else:
        actual_exon = first
        try:
            attrs = dict_exonCoords[ensg][actual_exon]
        except KeyError:
            chrom = 'unkonwn'
            strand = 'unknown'
            start_coord = 'unknown'
        else:
            chrom = attrs[0]
            strand = attrs[1]
            if strand == '+':
                start_coord = attrs[3]  # end
            else:
                start_coord = attrs[2]  # start
    
    # figure out end_coord
    if '_' in second:
        actual_exon,trailing = second.split('_')
        try:
            attrs = dict_exonCoords[ensg][actual_exon]
        except KeyError:
            if 'U' in actual_exon:
                end_coord = trailing
            elif 'ENSG' in actual_exon:
                end_coord = trailing
            else:
                end_coord = 'unknown'      
        else:
            end_coord = trailing
    else:
        actual_exon = second
        try:
            attrs = dict_exonCoords[ensg][actual_exon]
        except KeyError:
            if 'ENSG' in actual_exon:
                ensg_second, actual_exon_second = actual_exon.split(':')
                try:
                    attrs = dict_exonCoords[ensg_second][actual_exon_second]
                except KeyError:
                    end_coord = 'unknown'
                else:
                    if strand == '+':
                        end_coord = attrs[2]  # start
                    else:
                        end_coord = attrs[3]  # end
            else:
                end_coord = 'unknown'
        else:
            if strand == '+':
                end_coord = attrs[2]  # start
            else:
                end_coord = attrs[3]  # end
    
    # assemble
    if strand == '+':
        assemble = '{}:{}-{}({})'.format(chrom,start_coord,end_coord,strand)
    else:
        assemble = '{}:{}-{}({})'.format(chrom,end_coord,start_coord,strand)

    return assemble





