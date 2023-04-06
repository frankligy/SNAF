import pandas as pd
import numpy as np
import os,sys
import re
from tqdm import tqdm
import argparse
import multiprocessing as mp


def build_gene_model(df):
    # column names are: chr, source, type, start, end, score, strand, phase, attrs
    transcript_id = -1
    position_footprint = {}
    position_ense = {}
    gene_start, gene_end, strand, ensg, chromsome = None,None,None,None,None
    pat = re.compile(r'exon_id "(ENSE\d+)"')
    # phase1: build footprint
    for row in df.itertuples(index=False):
        if row.type == 'gene':
            gene_start, gene_end, strand, ensg, chromosome = row.start, row.end, row.strand, row.gene, row.chr
        elif row.type == 'transcript':
            transcript_id += 1
            continue
        elif row.type == 'exon':
            attrs = row.attrs
            ense = re.search(pat,attrs).group(1)
            for p in range(row.start,row.end+1,1):
                # put into position_footprint
                position_footprint.setdefault(p,[]).append(transcript_id)
                # put into position_ense
                position_ense.setdefault(p,[]).append(ense)
    position_ense = {p:list(set(ense_lis)) for p,ense_lis in position_ense.items()}
    # phase2: moving along the footprint
    if strand == '+':
        string_stream = ''
        running_profile = position_footprint[gene_start]
        block_index = 1
        intron_block_index = 0
        segment_index = 1
        anchor_position = gene_start
        p = gene_start
        while p <= gene_end:
            try:
                position_footprint[p]
            except KeyError:
                # write
                subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
                associated_ense = '|'.join(position_ense[p-1])
                string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg,chromosome,strand,subexon_identifier,anchor_position,p-1,associated_ense)  # because current p is the first base in intron, shouldn't be included in previous exon
                # enumerate the whole intron region
                pi = p
                while True:
                    try:
                        position_footprint[pi]
                    except KeyError:
                        pi += 1
                    else:
                        # write the intron
                        intron_block_index += 1
                        subexon_identifier = 'I'+str(intron_block_index)+'.'+str(1)
                        string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg,chromosome,strand,subexon_identifier,p,pi-1,'')  # because the current pi is the first base in next exon, shouldn't be included in the intron
                        # update
                        p = pi
                        running_profile = position_footprint[p]
                        block_index += 1
                        segment_index = 1
                        anchor_position = p
                        break
            else:
                if position_footprint[p] != running_profile:   
                    # write          
                    subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
                    associated_ense = '|'.join(position_ense[p-1])  # because the current p is the first base in next segment, we are writing the previous segment and its associated ense
                    string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg,chromosome,strand,subexon_identifier,anchor_position,p-1,associated_ense)
                    # update
                    running_profile = position_footprint[p]
                    segment_index += 1
                    anchor_position = p
                    p += 1
                else:
                    p += 1
        # write the last exon segment
        subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
        associated_ense = '|'.join(position_ense[p-1])
        string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg,chromosome,strand,subexon_identifier,anchor_position,p-1,associated_ense)

    else:   # negative string, so backtrack
        string_stream = ''
        running_profile = position_footprint[gene_end]
        block_index = 1
        intron_block_index = 0
        segment_index = 1
        anchor_position = gene_end
        p = gene_end
        while p >= gene_start:
            try:
                position_footprint[p]
            except KeyError:
                # write
                subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
                associated_ense = '|'.join(position_ense[p+1])
                string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg,chromosome,strand,subexon_identifier,p+1,anchor_position,associated_ense)  
                # enumerate the whole intron region
                pi = p
                while True:
                    try:
                        position_footprint[pi]
                    except KeyError:
                        pi -= 1
                    else:
                        # write the intron
                        intron_block_index += 1
                        subexon_identifier = 'I'+str(intron_block_index)+'.'+str(1)
                        string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg,chromosome,strand,subexon_identifier,pi+1,p,'')  
                        # update
                        p = pi
                        running_profile = position_footprint[p]
                        block_index += 1
                        segment_index = 1
                        anchor_position = p
                        break
            else:
                if position_footprint[p] != running_profile:   
                    # write          
                    subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
                    associated_ense = '|'.join(position_ense[p+1])  
                    string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg,chromosome,strand,subexon_identifier,p+1,anchor_position,associated_ense)
                    # update
                    running_profile = position_footprint[p]
                    segment_index += 1
                    anchor_position = p
                    p -= 1
                else:
                    p -= 1
        # write the last exon segment
        subexon_identifier = 'E'+str(block_index)+'.'+str(segment_index)
        associated_ense = '|'.join(position_ense[p+1])
        string_stream += '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(ensg,chromosome,strand,subexon_identifier,p+1,anchor_position,associated_ense)        
    return string_stream

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

def process_single_core(chunk):
    string_stream = ''
    for sub_df in tqdm(chunk,total=len(chunk)):
        sub_string_stream = build_gene_model(sub_df)
        string_stream += sub_string_stream
    return string_stream


def main(args):
    gtf = args.gtf
    gene = args.gene
    outdir = args.outdir
    df = pd.read_csv(gtf,sep='\t',header=None,skiprows=[0,1,2,3,4])
    df.columns = ['chr','source','type','start','end','score','strand','phase','attrs']
    df = df.loc[df['type'].isin(['gene','transcript','exon']),:]
    pat = re.compile(r'gene_id "(ENSG\d+)"')
    df['gene'] = [re.search(pat,item).group(1) for item in df['attrs']]
    if gene != 'all':
        df = df.loc[df['gene']==gene,:]
        string_stream = build_gene_model(df)
    else:
        string_stream = ''
        sub_df_list = list(list(zip(*list(df.groupby(by='gene'))))[1])
        n_cores = mp.cpu_count()
        chunks = split_array_to_chunks(sub_df_list,n_cores)
        pool = mp.Pool(processes=n_cores)
        print('spawn {} subprocesses'.format(n_cores))
        r = [pool.apply_async(func=process_single_core,args=(chunk,)) for chunk in chunks]
        pool.close()
        pool.join()
        for collect in r:
            sub_string_stream = collect.get()
            string_stream += sub_string_stream

    with open(os.path.join(outdir,'gene_model_{}.txt'.format(gene)),'w') as f:
        f.write(string_stream)
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build Gene Model from GTF')
    parser.add_argument('--gtf',type=str,default=None,help='the path to the gtf file')
    parser.add_argument('--gene',type=str,default=None,help='either all or stable ENSG ID')
    parser.add_argument('--outdir',type=str,default=None,help='output dir for the gene model txt file')
    args = parser.parse_args()
    main(args)


   

            








            
        







