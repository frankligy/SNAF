#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

import pandas as pd
import numpy as np
import os
import sys
import subprocess
import anndata as ad
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from sklearn.preprocessing import MinMaxScaler
import umap
import matplotlib.patches as mpatches
import subprocess
import pickle
import bisect
import re
from tqdm import tqdm



def big_dict_build(db_path):
    valid_chr = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13',
       'chr14','chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX','chrY']
    # we need a function which read in 'ENSG:E4.3-E5.1', output the coordinate and strand
    exons = pd.read_csv(db_path,sep='\t')
    exons = exons.loc[:,['ensembl_gene_id','exon_id','chromosome','strand','exon-region-start(s)','exon-region-stop(s)']]
    exons = exons.loc[exons['chromosome'].isin(valid_chr),:]
    exons.set_index(keys=pd.Index(np.arange(exons.shape[0])),inplace=True)

    # add one column called exon
    exon_col = [item.split('.')[0] for item in exons['exon_id']]
    exons.insert(loc=1,column='exon',value=exon_col)

    # start to build a query dict
    big_dict = {}  # { ENSG:(+,chr17,per_gene_dict), }
    for gene,subdf in tqdm(exons.groupby(by='ensembl_gene_id')):
        per_gene_dict = {}  # {'E1':coord, 'E2':coord}, coord is (start,end), it is all based on the order of forward string
        strand = subdf.iloc[0]['strand']  # which strand this gene lies on
        chromosome = subdf.iloc[0]['chromosome']  # which chromosome this gene lies on
        for exon,subsubdf in subdf.groupby(by='exon'):   
            subsubdf.sort_values(by='exon-region-start(s)',inplace=True)
            start = subsubdf.iloc[0]['exon-region-start(s)']   # strand doesn't matter here
            end = subsubdf.iloc[-1]['exon-region-stop(s)']
            coord = (start,end)
            per_gene_dict[exon] = coord
        big_dict[gene] = (strand,chromosome,per_gene_dict)
    return big_dict

def find_last_exon(exon):
    if exon.startswith('E'): # it is an exon
        try:
            assert int(exon[1:])-1 >= 1
        except AssertionError:   # means this exon, E1, is already the first one
            last = exon
        else:
            last = 'I' + str(int(exon[1:])-1)
    elif exon.startswith('I'): # it is an intron
        last = 'E' + exon[1:]
    return last

def find_next_exon(exon):
    '''
    One caveat for this function, E5 might be the last exon in the transcript
    '''
    if exon.startswith('E'):  # it is an exon
        next_ = 'I' + exon[1:]
    elif exon.startswith('I'):  # it is an intron
        next_ = 'E' + str(int(exon[1:]) + 1)
    return next_

def return_coords(big_dict,gene,query,start,end):
    '''
    big_dict stores the coord for every exon, not subexon
    gene is the ENSG id
    query is the query exon
    start is the exon that precedes the query exon, getting from find_last_exon function
    end is the exon that follows the query exon, getting from find_next_exon function
    '''
    global error   # a exon that is not present in db, 
    strand = big_dict[gene][0]
    chromosome = big_dict[gene][1]
    if strand == '+':
        try:
            final_1 = big_dict[gene][2][start][0]    # left-most
        except:
            error += 1
            return [1]
        try:
            final_2 = big_dict[gene][2][end][1]   # right-most
        except KeyError:   # E5 is the last exon, so if I5, it will throw a KeyError
            try:
                final_2 = big_dict[gene][2][query][1]   
            except:
                error += 1
                return [1]
    else:
        try:
            final_1 = big_dict[gene][2][end][0]
        except KeyError:
            try:
                final_1 = big_dict[gene][2][query][0]
            except:
                error += 1
                return [1]
        try:
            final_2 = big_dict[gene][2][start][1]
        except:
            error += 1
            return [1]
    '''return final1 always smaller than final2, because we are talking about the absolute position in positive strand'''
    return chromosome,strand,final_1,final_2

def build_peak_map(dir_path,filename):
    # assuming filename is like AARS.bed
    whole_path = os.path.join(dir_path,filename)
    IDR_peak = pd.read_csv(whole_path,sep='\t',header=None)
    IDR_peak = IDR_peak.sort_values(by=11,ascending=False)   # global IDR -log10
    IDR_peak.set_index(keys=np.arange(IDR_peak.shape[0]),inplace=True)
    IDR_peak['percentile'] = IDR_peak.index.values / IDR_peak.shape[0]
    peak_map = {}  
    peak_sig = {}
    for chromosome, subdf in IDR_peak.groupby(by=0):
        valid_chrom = ['chr{}'.format(i) for i in range(1,23)] + ['chr{}'.format(i) for i in ['X','Y']]
        if chromosome in valid_chrom:
            peak_map[chromosome] = {'+':[],'-':[]}
            peak_sig[chromosome] = {'+':[],'-':[]}
            for strand,sub2df in subdf.groupby(by=5):
                sub2df = sub2df.sort_values(by=1)
                start_pos = sub2df[1].values
                end_pos = sub2df[2].values
                sig_value = 1 - sub2df['percentile'].values
                for s,e,v in zip(start_pos,end_pos,sig_value):
                    peak_map[chromosome][strand].extend([s,e])
                    peak_sig[chromosome][strand].append(v)
                sorted_coord = sorted(peak_map[chromosome][strand])
                try:
                    assert peak_map[chromosome][strand] == sorted_coord
                except AssertionError:
                    raise Exception(chromsome,strand,'order incorrect!!') 
    total_peak_map[filename.split('.')[0]] = peak_map
    total_peak_sig[filename.split('.')[0]] = peak_sig

def disentangle_uid(uid,add_stable_id=False,add_symbol=False):
    fg,bg = uid.split('|')
    gene_symbol = fg.split(':')[0]
    stable_id = fg.split(':')[1]
    exons = ':'.join(fg.split(':')[2:])
    output = exons
    if add_stable_id:
        output = stable_id + ':' + output
    if add_symbol:
        output = gene_symbol + ':' + output
    return output


# start main program
ee_path = '/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/AltResults/AlternativeOutput/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt'
event_and_alt = pd.read_csv(ee_path,sep='\t').loc[:,['UID','AltExons']].set_index('UID').squeeze()  # series
event_and_alt.index = [disentangle_uid(item,add_symbol=True) for item in event_and_alt.index]
event_and_alt = event_and_alt.loc[np.logical_not(event_and_alt.index.duplicated())]


# build or load big_dict
# build
# big_dict = big_dict_build('../Hs_Ensembl_exons_91.txt')
# with open('../big_dict.p','wb') as f:
#     pickle.dump(big_dict,f,protocol=4)
# load
with open('../big_dict.p','rb') as f:
    big_dict = pickle.load(f)
    # { ENSG:(+,chr17,per_gene_dict), }
    # per_gene_dict {'E1':coord, 'E2':coord}, coord is (start,end), it is all based on the order of forward string


bed_folder_path = '/data/salomonis2/NCI-R01/Signatures/encode_eclip_bed/raw/HepG2/idr_bed'
ori_dir = os.getcwd()  
os.chdir(bed_folder_path)
all_beds = subprocess.run("for file in *.bed; do echo $file; done",shell=True,stdout=subprocess.PIPE,universal_newlines=True).stdout.split('\n')[:-1]
os.chdir(ori_dir)

total_peak_map = {} # {U2AF2:{chr1:{'+':[1,4,7]}},PTBP1:{chr1:[4,5,7]}}   each two coords corresponds to a sig value
total_peak_sig = {}  # {U2AF2:{chr1:{'+':[0.5,0.6]}}}}
for bed in tqdm(all_beds):
    build_peak_map(bed_folder_path,bed)

# start to build priors
gs = [] # [(sf,event,percentile)]
error = 0
for key,value in tqdm(total_peak_map.items()):
    sig_value = total_peak_sig[key]
    for event,exons in zip(event_and_alt.index,event_and_alt.values):
        exonlist = exons.split('|')
        for exon in exonlist:
            pattern = re.compile(r'^ENSG\d+:[EI]\d+\.\d+$')  # don't consider complicated ones
            if re.search(pattern,exon):
                gene = exon.split(':')[0]
                query = exon.split(':')[1].split('.')[0]  # exon not sub-exon, so it should be E4 instead of E4.2
                start = find_last_exon(query)
                end = find_next_exon(query)
                inspect = return_coords(big_dict, gene, query, start, end)  # inspect (chrom,strand,final_1,final_2)
                if len(inspect) == 1: # error
                    continue
                chr_ = inspect[0]
                strand = inspect[1]
                try:
                    index1 = bisect.bisect(value[chr_][strand],inspect[2])
                    index2 = bisect.bisect(value[chr_][strand],inspect[3])
                except KeyError: # if peak_map doesn't have chr20, it will have keyError, just continue to next one
                    continue
                else:
                    if index2 - index1 >= 1:  # at least one peak falls into the interval
                        max_associated_sig = 0
                        if index1 % 2 == 0: # meaning the left end falls into inter-peak
                            left_sig_index = index1 / 2
                            if index2 % 2 == 0: # meaning the right end falls into inter-peak
                                right_sig_index = index2 / 2 - 1
                            else:
                                right_sig_index = (index2-1)/2
                        else:
                            left_sig_index = (index1-1)/2
                            if index2 % 2 == 0:
                                right_sig_index = index2 / 2 - 1
                            else:
                                right_sig_index = (index2-1)/2
                        for sig_index in range(int(left_sig_index),int(right_sig_index)+1,1):
                            associated_sig = sig_value[chr_][strand][sig_index]
                            if associated_sig > max_associated_sig:
                                max_associated_sig = associated_sig               
                        gs.append((key,event,max_associated_sig))



 
# make the gs a dataframe, aka, prior matrix
prior_binding = pd.DataFrame()
tmp = list(zip(*gs))
prior_binding['sf'] = tmp[0]
prior_binding['event'] = tmp[1]
prior_binding['value'] = tmp[2]
prior_binding= prior_binding.groupby(by=['event','sf'])['value'].mean().unstack(fill_value=0)
prior_binding.to_csv('/data/salomonis-archive/FASTQs/NCI-R01/ENCODE_KD_fastq_frank/altanalyze_output/prior/binding_prior/HepG2_prior.txt',sep='\t')




