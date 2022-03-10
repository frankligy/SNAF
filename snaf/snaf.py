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




# configuration
def snaf_configuration(exon_table,fasta,software_path_arg=None,binding_method_arg=None):
    global dict_exonCoords
    global dict_fa
    global software_path
    global binding_method
    dict_exonCoords = exonCoords_to_dict(exon_table)
    dict_fa = fasta_to_dict(fasta)
    software_path = software_path_arg
    binding_method = binding_method_arg



class JunctionCountMatrixQuery():
    def __init__(self,junction_count_matrix,cores=None,add_control=None):
        self.junction_count_matrix = junction_count_matrix
        if cores is None:
            cores = mp.cpu_count()
        self.cores = cores
        self.get_neojunctions(add_control=add_control)
        

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
    
    def get_neojunctions(self,add_control):
        self.valid, self.invalid, self.cond_df = multiple_crude_sifting(self.junction_count_matrix,add_control)
        self.subset = self.junction_count_matrix.loc[self.valid,:]
        self.cond_subset_df = self.cond_df.loc[self.valid,:]

    def get_neojunction_info(self,event):  # more for B cell antigen, since we design it in a way that B cell pipeline completely rely on T pipeline for GTEx check
        ed = self.subset.loc[event,:].to_dict()
        freq = np.count_nonzero(np.array(list(ed.values())))/len(ed)
        return ed,freq


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
    def each_chunk_func(input_,kind,hlas=None,sub_cond_df=None,binding_method=None):
        if kind == 1:
            nj_list = []
            for uid in input_.index:
                nj = NeoJunction(uid=uid,count=50,check_gtex=False)
                nj.detect_type()
                nj.retrieve_junction_seq()
                nj.in_silico_translation()    
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

    def run(self,hlas,outdir='.',name='after_prediction.p'):
        self.parallelize_run(kind=1)
        self.parallelize_run(kind=3,hlas=hlas)
        self.serialize(outdir=outdir,name=name)

    @staticmethod
    def generate_results(self,path,outdir):
        jcmq = snaf.JunctionCountMatrixQuery.deserialize(name=path)
        snaf.downstream.stage0_compatible_results(jcmq,outdir=outdir)
        for stage in [3,2,1]:
            jcmq.show_neoantigen_burden(outdir=outdir,name='burden_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False)
            jcmq.show_neoantigen_frequency(outdir=outdir,name='frequency_stage{}.txt'.format(stage),stage=stage,verbosity=1,contain_uid=False,plot=True,plot_name='frequency_stage{}.pdf'.format(stage))
            jcmq.show_neoantigen_frequency(outdir=outdir,name='frequency_stage{}_verbosity1_uid.txt'.format(stage),stage=stage,verbosity=1,contain_uid=True,plot=False)

    def parallelize_run(self,kind,hlas=None):
        pool = mp.Pool(processes=self.cores)
        if kind == 1 or kind == 2:
            sub_dfs = JunctionCountMatrixQuery.split_df_to_chunks(self.subset,self.cores)
            r = [pool.apply_async(func=JunctionCountMatrixQuery.each_chunk_func,args=(sub_df,kind,)) for sub_df in sub_dfs]
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
            r = [pool.apply_async(func=JunctionCountMatrixQuery.each_chunk_func,args=(sub_array,kind,hlas,sub_cond_df,binding_method,)) for sub_array,sub_cond_df in zip(sub_arrays,sub_cond_dfs)]
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

    def visualize(self,uid,sample,outdir,tumor=False):
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        if sample is not None:
            row_index = self.subset.index.tolist().index(uid)
            col_index = self.subset.columns.tolist().index(sample)
            results = self.results[0]
            hlas = self.results[1]
            nj = deepcopy(results[row_index])
            nj.enhanced_peptides = nj.enhanced_peptides.filter_based_on_hla(selected_hla=hlas[col_index])
            nj.visualize(outdir,'{}_{}.pdf'.format(uid.replace(':','_'),sample))
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




    def show_neoantigen_burden(self,outdir,name,stage,verbosity,contain_uid):
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        sub_arrays = JunctionCountMatrixQuery.split_array_to_chunks(self.results[0],self.cores)
        sub_conds = JunctionCountMatrixQuery.split_df_to_chunks(self.cond_subset_df,self.cores)
        hlas = self.results[1]
        pool = mp.Pool(processes=self.cores)

        r = [pool.apply_async(func=JunctionCountMatrixQuery.show_neoantigen_burden_single_run,args=(sub_array,sub_cond,hlas,stage,verbosity,contain_uid,)) for sub_array,sub_cond in zip(sub_arrays,sub_conds)]
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

    def show_neoantigen_frequency(self,outdir,name,stage,verbosity,contain_uid,plot,plot_name=None,yscale='linear'):
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        sub_arrays = JunctionCountMatrixQuery.split_array_to_chunks(self.results[0],self.cores)
        sub_conds = JunctionCountMatrixQuery.split_df_to_chunks(self.cond_subset_df,self.cores)
        hlas = self.results[1]
        column_names = self.subset.columns
        pool = mp.Pool(processes=self.cores)

        r = [pool.apply_async(func=JunctionCountMatrixQuery.show_neoantigen_frequency_single_run,args=(sub_array,sub_cond,hlas,column_names,stage,verbosity,contain_uid,)) for sub_array,sub_cond in zip(sub_arrays,sub_conds)]
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
            sns.histplot(df['n_sample'].values,binwidth=1,kde=True,ax=ax)
            ax.set_yscale(yscale)
            plt.savefig(os.path.join(outdir,'x_occurence_{}'.format(plot_name)),bbox_inches='tight')
            plt.close()




    def show_neoantigen_as_fasta(self,outdir,name,stage,verbosity,contain_uid,sample=None):
        if not os.path.exists(outdir):
            os.mkdir(outdir)
            
        sub_arrays = JunctionCountMatrixQuery.split_array_to_chunks(self.results[0],self.cores)
        sub_conds = JunctionCountMatrixQuery.split_df_to_chunks(self.cond_subset_df,self.cores)
        hlas = self.results[1]
        column_names = self.subset.columns
        pool = mp.Pool(processes=self.cores)

        r = [pool.apply_async(func=JunctionCountMatrixQuery.show_neoantigen_frequency_single_run,args=(sub_array,sub_cond,hlas,column_names,stage,verbosity,contain_uid,)) for sub_array,sub_cond in zip(sub_arrays,sub_conds)]
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
    def show_neoantigen_burden_single_run(sub_array,sub_cond,hlas,stage,verbosity,contain_uid):
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
                    nj_copy.derive_candidates(stage=stage,verbosity=verbosity,contain_uid=contain_uid)
                    nj_burden.append(len(nj_copy.candidates))
            nj_burden_2d.append(nj_burden)  
        df = pd.DataFrame(data=nj_burden_2d)
        return df

    # let's define the single function to run in each subprocess
    @staticmethod
    def show_neoantigen_frequency_single_run(sub_array,sub_cond,hlas,column_names,stage,verbosity,contain_uid):
        dic = {}
        for nj,index in zip(sub_array,range(sub_cond.shape[0])):
            cond_row = sub_cond.iloc[index].values
            for hla,column_name,cond in zip(hlas,column_names,cond_row):
                if nj is not None and cond:
                    nj_copy = deepcopy(nj)
                    nj_copy.enhanced_peptides = nj_copy.enhanced_peptides.filter_based_on_hla(selected_hla=hla)
                    nj_copy.derive_candidates(stage=stage,verbosity=verbosity,contain_uid=contain_uid)
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
        if kind == 2:  # hla-reduced enhanced peptide
            self.mers = peptides
            self.info = hlas
        else:
            self.mers = []
            self.info = []
            for k,v in peptides.items():
                self.mers.append(k)
                phlas = {}  # {pep1:{origin:(extra,n_from_first),hla1:{},hla2:{}}}
                for pep in v:
                    if kind == 0:  # each pep is (pep,extra,n_from_first)
                        pairs = {}
                        pairs['origin'] = (pep[1],pep[2])
                        for hla in hlas:
                            pairs[hla] = {}
                        phlas[pep[0]] = pairs
                    elif kind == 1:   # each pep is (pep,extra,n_from_first,hla)
                        try:
                            phlas[pep[0]]['origin'] = (pep[1],pep[2])
                        except KeyError:
                            phlas[pep[0]] = {}
                            phlas[pep[0]]['origin'] = (pep[1],pep[2])
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
                    extra,n_from_first = hla_complex['origin']
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
                extra,n_from_first = hla_complex['origin']
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
                        peptides[k].append((pep,extra,n_from_first,hla))   # if not reinstantiate, this format is called reduced form
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
            seq1 = subexon_tran(subexon1,ensid,'site1')
            seq2 = subexon_tran(subexon2,ensid,'site2')
            junction = ','.join([seq1,seq2])
            self.junction = junction
        else:
            self.junction = '$' * 10   # indicating invalid uid

    def in_silico_translation(self,ks=[9,10]):
        peptides = {k:[] for k in ks}
        if '$' not in self.junction and '*' not in self.junction and '#' not in self.junction:
            first,second = self.junction.split(',')
            for phase in [0,1,2]:  # tranlation starts with index "phase"
                de_facto_first = first[phase:]
                pep_dict = get_peptides(de_facto_first,second,ks)
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


    def derive_candidates(self,stage,verbosity,contain_uid=True):
        if stage == 1: # all translated peptides
            self.candidates = self.enhanced_peptides.simplify_to_list(verbosity=verbosity)
        elif stage == 2: # all bound peptides
            self.candidates = self.enhanced_peptides.filter_based_on_criterion([('netMHCpan_el',0,'<=',2),]).simplify_to_list(verbosity=verbosity)
        elif stage == 3: # immunogenic peptides
            self.candidates = self.enhanced_peptides.filter_based_on_criterion([('netMHCpan_el',0,'<=',2),('deepimmuno_immunogenicity',1,'==','True'),]).simplify_to_list(verbosity=verbosity)
        if contain_uid:
            new = []
            for item in self.candidates:
                tmp = [str(i) for i in item]
                tmp.append(self.uid)
                new.append(','.join(tmp))
            self.candidates = new





    def visualize(self,outdir,name):
        reduced = self.enhanced_peptides.filter_based_on_criterion([('netMHCpan_el',0,'<=',2),('deepimmuno_immunogenicity',1,'==','True'),],False)
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
                aa, extra, n_from_first, hla = to_draw[i]
                first,second = self.junction.split(',')
                dna_first = extra + n_from_first * 3
                dna_second = -extra + (len(aa)-n_from_first) * 3
                binding_score = self.enhanced_peptides[len(aa)][aa][hla]['netMHCpan_el'][0]
                immunogenicity_score = self.enhanced_peptides[len(aa)][aa][hla]['deepimmuno_immunogenicity'][0]   
                # draw     
                ax = show_candicates(ax,aa,extra,n_from_first,hla,first,second,dna_first,dna_second,binding_score,immunogenicity_score) 
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





def get_peptides(de_facto_first,second,ks):
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
                        peptides[k].append((pep,extra,n_from_first))
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

def subexon_tran(subexon,EnsID,flag):  # flag either site1 or site2
    '''
    1. subexon can take multiple forms depending on the event type
    E1.2 or I3.4
    E6.2_67996641 or I40.1_13076665, also depending on whether they are subexon1 or subexon2
    ENSG00000213934:E3.1 or ENSG00000213934:E2.1_473843893894
    U0.1_49689185
    2. everything with trailing suffix will depend on the subexon1 or subexon2, but sometimes, it is fixed (trans-splicing can only be in subexon2)
    3. to be clear, the exon_seq returned is always 5'-3' sequence, not forward anymore.
    '''
    try:   # E1.2 or I3.4
        attrs = dict_exonCoords[EnsID][subexon]  # [chr,strand,start,end,suffer] 
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


def add_coord_frequency_table(df,remove_quote=True):
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
                attrs = dict_exonCoords[ensg_second][actual_exon_second]
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





