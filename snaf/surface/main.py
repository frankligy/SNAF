import pandas as pd
import numpy as np
import os,sys
import subprocess
import multiprocessing as mp
import re
import pickle
from datetime import datetime,date
from tqdm import tqdm
from .data_io import *
from .orf_finder import *
from .orf_check import *
from .alignment import *



def initialize(transcript_db,exon_table,fasta,membrane_db,biotype_db,membrane_fasta_db):
    global df_exonlist
    global dict_exonCoords
    global dict_fa
    global dict_biotype
    global df_membrane_proteins
    global dict_uni_fa
    print('{} {} starting surface antigen initialization'.format(date.today(),datetime.now().strftime('%H:%M:%S')))
    df_exonlist = pd.read_csv(transcript_db,sep='\t',header=None,names=['EnsGID','EnsTID','EnsPID','Exons'])  # index is number
    dict_exonCoords = exonCoords_to_dict(exon_table) 
    dict_fa = fasta_to_dict(fasta)
    dict_biotype = biotype(pd.read_csv(biotype_db,sep='\t'))  # index is number
    df_membrane_proteins = pd.read_csv(membrane_db,sep='\t',index_col=0)
    dict_uni_fa = read_uniprot_seq(membrane_fasta_db)
    print('{} {} finished surface antigen initialization'.format(date.today(),datetime.now().strftime('%H:%M:%S')))

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

def single_run(uids,n_stride,tmhmm,software_path):
    results = []
    for uid in tqdm(uids):
        sa = SurfaceAntigen(uid,False)
        sa.detect_type()
        sa.retrieve_junction_seq()
        sa.recovery_full_length_protein()
        sa.find_orf()
        sa.orf_check(n_stride=n_stride)
        sa.align_uniprot(tmhmm=tmhmm,software_path=software_path)
        results.append(sa)
    return results

def batch_run(uid_list,cores,n_stride,tmhmm=False,software_path=None,outdir='.',name=None):
    sub_uid_lists = split_array_to_chunks(uid_list,cores=cores)
    pool = pool = mp.Pool(processes=cores)
    r = [pool.apply_async(func=single_run,args=(sub_uid_list,n_stride,tmhmm,software_path,)) for sub_uid_list in sub_uid_lists]
    pool.close()
    pool.join()
    results = []
    for collect in r:
        result = collect.get()
        resutls.extend(result)
    if name is None:
        name = 'batch_run_surface_antigen.p'
    with open(os.path.join(outdir,name),'wb') as f:
        pickle.dump(results,f)


def get_all_transcripts(ensgid,outdir='.'):
    df_certain = df_exonlist.loc[df_exonlist['EnsGID']==ensgid,:]
    df_certain.to_csv(os.path.join(outdir,'{}_all_transcripts.txt'.format(ensgid)),sep='\t',index=None)


def get_existing_isoforms(ensgid,outdir='.'):
    with open(os.path.join(outdir,'{}_existing_isoforms.fasta'.format(ensgid)),'w') as f:
        for k,v in dict_uni_fa[ensgid].items():
            f.write('>{}\n{}\n'.format(k,v))

def filter_to_membrane_protein(lis):
    filtered_lis = []
    all_membrane = set(dict_uni_fa.keys())
    for uid in lis:
        ensgid = uid.split(':')[0]
        if ensgid in all_membrane:
            filtered_lis.append(uid)
    return filtered_lis



class SurfaceAntigen(object):

    def __init__(self,uid,check_overlap=True):
        self.uid = uid
        if check_overlap:
            if not self.is_membrane_protein():
                raise Exception('This event will not encode a surface protein')

    def is_membrane_protein(self):
        ensgid = self.uid.split(':')[0]
        if ensgid in set(df_membrane_proteins['Ens'].tolist()):
            return True
        else:
            return False

    def __str__(self):
        print_str = 'uid:{}\n'.format(self.uid)
        try:
            print_event_type = self.event_type
        except AttributeError:
            print_event_type = None
        print_str += 'event type:{}\n'.format(print_event_type)
        try:
            print_junction = self.junction[:5] + '...' + self.junction[-5:]
        except AttributeError:
            print_junction = None
        print_str += 'Junction:{}\n'.format(print_junction)   
        try:
            if self.full_length == ['unrecoverable']:
                print_full_length = (len(self.full_length),self.full_length)
            else:
                print_full_length = (len(self.full_length),[self.full_length.index(item) for item in self.full_length if item != ''])
        except AttributeError:
            print_full_length = (None,None)
        print_str += 'Full length transcripts: length {}, indices {}\n'.format(print_full_length[0],print_full_length[1])
        try:
            if self.orft == ['unrecoverable']:
                print_orft = (len(self.orft),self.orft)
            else:
                print_orft = (len(self.orft),[self.orft.index(item) for item in self.orft if item != ''])
        except AttributeError:
            print_orft = (None,None)
        print_str += 'ORF transcripts: length {}, indices {}\n'.format(print_orft[0],print_orft[1])
        try:
            if self.orfp == ['unrecoverable']:
                print_orfp = (len(self.orfp),self.orfp)
            else:
                print_orfp = (len(self.orfp),[self.orfp.index(item) for item in self.orfp if item != ''])
        except AttributeError:
            print_orfp = (None,None)
        print_str += 'ORF peptides: length {}, indices {}\n'.format(print_orfp[0],print_orfp[1])
        try:
            if self.nmd == ['unrecoverable']:
                print_nmd = (len(self.nmd),self.nmd)
            else:
                print_nmd = (len(self.nmd),[item for item in self.nmd if item != ''])
        except AttributeError:
            print_nmd =(None,None)  
        print_str += 'NMD check: length {}, indices {}\n'.format(print_nmd[0],print_nmd[1])      
        try:
            if self.translatability == ['unrecoverable']:
                print_translatability = (len(self.translatability),self.translatability)
            else:
                print_translatability = (len(self.translatability),[item for item in self.translatability if item != ''])
        except AttributeError:
            print_translatability = (None,None)
        print_str += 'tranlatability check: length {}, indices {}\n'.format(print_translatability[0],print_translatability[1])    
        try:
            if self.alignment == ['unrecoverable']:
                print_alignment = (len(self.alignment),self.alignment)
            else:
                print_alignment = (len(self.alignment),[item for item in self.alignment if item != ''])
        except AttributeError:
            print_alignment = (None,None)
        print_str += 'Alignment: length {}, indices {}\n'.format(print_alignment[0],print_alignment[1])   
        return print_str



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


    def recovery_full_length_protein(self):
        if '$' not in self.junction and '*' not in self.junction and '#' not in self.junction:
            ensgid = self.uid.split(':')[0]
            exons = ':'.join(self.uid.split(':')[1:])
            if self.event_type == 'ordinary':
                full_transcript_store = first_round_match_novel(ensgid,exons)
            elif self.event_type == 'alt3' or self.event_type == 'alt5' or self.event_type == 'alt3_alt5':
                full_transcript_store = second_round_match(ensgid,exons)
            elif self.event_type == 'intron_retention':
                full_transcript_store = third_round_match(ensgid,exons)
            elif self.event_type == 'trans_splicing':
                full_transcript_store = fourth_round_match(ensgid,exons)
            else:   # novel_exon and utr_event
                full_transcript_store = ['unrecoverable']
        else:
            full_transcript_store = ['unrecoverable']
        self.full_length = full_transcript_store

    def find_orf(self):
        orft_list = []
        orfp_list = []
        for sequence in self.full_length:
            if sequence != 'unrecoverble' and sequence != '':
                candidate_orfs = transcript2orf(sequence)
                max_orf = prioritize_orf(candidate_orfs)
                max_pep = orf2pep(max_orf)
                orft_list.append(max_orf)
                orfp_list.append(max_pep)
            else:
                orft_list.append(sequence)
                orfp_list.append(sequence)
        self.orft = orft_list
        self.orfp = orfp_list

    def orf_check(self,n_stride):
        set_global_env(df_exonlist,dict_exonCoords,dict_fa,dict_biotype)
        nmd_check_result = nmd_check(self.uid,self.full_length,self.orft,n_stride)
        translatability_check_result = translatability_check(self.uid,self.orft)
        self.nmd = nmd_check_result
        self.translatability = translatability_check_result

    def align_uniprot(self,tmhmm,software_path=None):
        results = alignment_to_uniprot(self.orfp,self.uid,dict_uni_fa,tmhmm,software_path)
        self.alignment = results
        if tmhmm:
            subprocess.run('rm -r ./TMHMM_*',shell=True)  

    def visualize(self,index,outdir='.',name=None,fragment=None):
        full_length = self.full_length[index]
        junction = self.junction.replace(',','')
        orft = self.orft[index]
        ensgid = self.uid.split(':')[0]
        exonlist = df_exonlist.loc[df_exonlist['EnsGID']==ensgid,:]['Exons'].iloc[index].split('|')
        if full_length == '' or orft == '':
            raise Exception('please select index that are not empty based on SurfaceAntigen summary')
        else:
            import matplotlib.pyplot as plt
            from matplotlib.patches import Rectangle,Patch
            l = len(full_length)
            start_j = full_length.index(junction)
            end_j = start_j + (len(junction)-1)
            start_o = full_length.index(orft)
            end_o = start_o + (len(orft)-1)
            fig,ax = plt.subplots()
            ax.set_xlim(-0.05,1.05)
            ax.set_ylim(-0.05,1.05)
            # draw full_length
            rect_full_length = Rectangle((0,0.8),1,0.1,linewidth=0.5,facecolor='g',edgecolor='k')
            ax.add_patch(rect_full_length)
            # draw orft
            rect_orft = Rectangle((start_o/l,0.6),(end_o-start_o)/l,0.1,linewidth=0.5,facecolor='orange',edgecolor='k')
            ax.add_patch(rect_orft)
            # draw junction
            rect_junction = Rectangle((start_j/l,0.4),(end_j-start_j)/l,0.1,linewidth=0.5,facecolor='r',edgecolor='k')
            ax.add_patch(rect_junction)
            # draw exonlist
            for i,exon in enumerate(exonlist):
                seq = get_exon_seq(exon,ensgid)
                try:
                    start_s = full_length.index(seq)
                except ValueError:   # say E9.2-E13.1 is the novel splicing event, skip the E11.1 in the list, so E11.1 won't align
                    continue 
                end_s = start_s + (len(seq)-1)
                if i % 2 == 0:
                    rect_seq = Rectangle((start_s/l,0.2),(end_s-start_s)/l,0.1,linewidth=0.1,facecolor='pink',edgecolor='k')
                else:
                    rect_seq = Rectangle((start_s/l,0.2),(end_s-start_s)/l,0.1,linewidth=0.1,facecolor='b',edgecolor='k')
                ax.add_patch(rect_seq)
                ax.text(x=(start_s + end_s)/2/l,y=0.3,s=exon,rotation=90,fontsize=2,va='bottom')
            # draw fragment
            if fragment is not None:
                start_f = full_length.index(fragment)
                end_f = start_f + (len(fragment)-1)
                rect_fragment = Rectangle((start_f/l,0.0),(end_f-start_f)/l,0.1,linewidth=0.5,facecolor='magenta',edgecolor='k')
                ax.add_patch(rect_fragment)
            # draw legend
            ax.legend(handles=[Patch(color=i) for i in ['g','orange','r','magenta']],labels=['transcript','ORF','junction','fragment'],
                      bbox_to_anchor=(1,1),loc='upper left',frameon=False)
            # draw title
            ax.set_title('{}_{}'.format(self.uid,index))
            if name is None:
                name = '{}_{}.pdf'.format(self.uid.replace(':','_'),index)
            plt.savefig(os.path.join(outdir,name),bbox_inches='tight')
            plt.close()

        


    


    
# standalone functions
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
        exon_seq = get_exon_seq(subexon,EnsID,check=False,attrs=attrs)
    except KeyError:
        if ':' in subexon: # ENSG00000213934:E3.1
            fusionGeneEnsID = subexon.split(':')[0] 
            fusionGeneExon = subexon.split(':')[1]        
            if '_' in fusionGeneExon:   # ENSG:E2.1_473843893894
                exon_seq = get_trailing_exon_seq(fusionGeneExon,fusionGeneEnsID,'site2',check=False)
            else:  # ENSG:E2.1
                try:
                    attrs = dict_exonCoords[fusionGeneEnsID][fusionGeneExon]
                except KeyError:
                    exon_seq = '*' * 10  # indicator for error on MultiPath-PSI itself
                else:
                    exon_seq = get_exon_seq(fusionGeneExon,fusionGeneEnsID,check=False,attrs=attrs)

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
                    exon_seq = get_trailing_exon_seq('_'.join([subexon,suffix]),EnsID,flag=flag,check=False)
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


def get_exon_seq(exon,ensgid,check=True,attrs=None):
    if check:
        try:
            attrs = dict_exonCoords[ensgid][exon]
        except KeyError:   # indicator for error on MultiPath-PSI itself
            frag = '*' * 10
        else:
            strand = attrs[1]
            suffer = attrs[4]
            if strand == '+' and suffer == 'False':
                frag = query_from_dict_fa(attrs[2],attrs[3],ensgid,attrs[1])
            elif strand == '+' and suffer == 'True':
                frag = query_from_dict_fa(int(attrs[2])-1,attrs[3],ensgid,attrs[1])
            elif strand == '-' and suffer == 'False':
                frag = query_from_dict_fa(attrs[2],attrs[3],ensgid,attrs[1])
            elif strand == '-' and suffer == 'True':
                frag = query_from_dict_fa(attrs[2],int(attrs[3])-1,ensgid,attrs[1])
    else:
        strand = attrs[1]
        suffer = attrs[4]
        if strand == '+' and suffer == 'False':
            frag = query_from_dict_fa(attrs[2],attrs[3],ensgid,attrs[1])
        elif strand == '+' and suffer == 'True':
            frag = query_from_dict_fa(int(attrs[2])-1,attrs[3],ensgid,attrs[1])
        elif strand == '-' and suffer == 'False':
            frag = query_from_dict_fa(attrs[2],attrs[3],ensgid,attrs[1])
        elif strand == '-' and suffer == 'True':
            frag = query_from_dict_fa(attrs[2],int(attrs[3])-1,ensgid,attrs[1])
    return frag

def get_trailing_exon_seq(exon,ensgid,flag,check=True):
    e,suffix = exon.split('_')
    if check:
        try:
            attrs = dict_exonCoords[ensgid][e]
        except KeyError:
            exon_seq = '*' * 10
        else:
            if flag == 'site2':
                if attrs[1] == '+':  
                    if attrs[4] == 'True': 
                        exon_seq = query_from_dict_fa(suffix,attrs[3],ensgid,attrs[1]) 
                    else:
                        exon_seq = query_from_dict_fa(suffix,attrs[3],ensgid,attrs[1]) 
                else:  
                    if attrs[4] == 'True': 
                        exon_seq = query_from_dict_fa(attrs[2],suffix,ensgid,attrs[1])
                    else:
                        exon_seq = query_from_dict_fa(attrs[2],suffix,ensgid,attrs[1])
            elif flag == 'site1':  
                if attrs[1] == '+': 
                    if attrs[4] == 'True':
                        exon_seq = query_from_dict_fa(int(attrs[2])-1,suffix,ensgid,attrs[1])
                    else:
                        exon_seq = query_from_dict_fa(attrs[2],suffix,ensgid,attrs[1])
                else:
                    if attrs[4] == 'True':
                        exon_seq = query_from_dict_fa(suffix,int(attrs[3])-1,ensgid,attrs[1])
    else:
        attrs = dict_exonCoords[ensgid][e]
        if flag == 'site2':
            if attrs[1] == '+':  
                if attrs[4] == 'True': 
                    exon_seq = query_from_dict_fa(suffix,attrs[3],ensgid,attrs[1]) 
                else:
                    exon_seq = query_from_dict_fa(suffix,attrs[3],ensgid,attrs[1]) 
            else:  
                if attrs[4] == 'True': 
                    exon_seq = query_from_dict_fa(attrs[2],suffix,ensgid,attrs[1])
                else:
                    exon_seq = query_from_dict_fa(attrs[2],suffix,ensgid,attrs[1])
        elif flag == 'site1':  
            if attrs[1] == '+': 
                if attrs[4] == 'True':
                    exon_seq = query_from_dict_fa(int(attrs[2])-1,suffix,ensgid,attrs[1])
                else:
                    exon_seq = query_from_dict_fa(attrs[2],suffix,ensgid,attrs[1])
            else:
                if attrs[4] == 'True':
                    exon_seq = query_from_dict_fa(suffix,int(attrs[3])-1,ensgid,attrs[1])
    return exon_seq




def first_round_match_existing(ensgid,exons):
    # exons is E4.1-E5.6
    # first round mainly solve ordinary event
    exons = exons.replace('-','|')
    full_transcript_store = []  # ['',full_transcript1_seq,...] 
    df_certain = df_exonlist.loc[df_exonlist['EnsGID']==ensgid,:]
    for item in list(df_certain['Exons']):
        full_transcript = ''
        pattern1 = re.compile(r'\|{}\|'.format(exons.replace('|','\|')))  # |E3.4|E4.5|
        pattern2 = re.compile(r'^{}\|'.format(exons.replace('|','\|')))    # E1.1|E4.5|  first one
        pattern3 = re.compile(r'\|{}$'.format(exons.replace('|','\|')))     # |E6.9|E7.1  last one 
        pattern4 = re.compile(r'^{}$'.format(exons.replace('|','\|')))     # E4.5|E5.6   just it
        if re.search(pattern1,item) or re.search(pattern2,item) or re.search(pattern3,item) or re.search(pattern4,item):
            exonlist = item.split('|')
            for exon in exonlist:
                frag = get_exon_seq(exon,ensgid)
                full_transcript += frag
            full_transcript_store.append(full_transcript)   
        else:
            full_transcript_store.append('')
    return full_transcript_store

def first_round_match_novel(ensgid,exons):
    # exons is E4.1-E5.6
    # match E4.1 and E5.6 separately
    e1,e2 = exons.split('-')
    full_transcript_store = []  # ['',full_transcript1_seq,...] 
    df_certain = df_exonlist.loc[df_exonlist['EnsGID']==ensgid,:]
    pattern1_1 = re.compile(r'{}\|'.format(e1))
    pattern1_2 = re.compile(r'{}$'.format(e1))
    pattern2_1 = re.compile(r'{}\|'.format(e2))
    pattern2_2 = re.compile(r'{}$'.format(e2))
    for item in list(df_certain['Exons']):
        match1 = re.search(pattern1_1,item) or re.search(pattern1_2,item)
        match2 = re.search(pattern2_1,item) or re.search(pattern2_2,item)
        if match1 and match2:
            full_transcript = ''
            exonlist = item.split('|')
            index1 = exonlist.index(e1)
            index2 = exonlist.index(e2)
            for i in range(index1):
                full_transcript += get_exon_seq(exonlist[i],ensgid)
            full_transcript += get_exon_seq(exonlist[index1],ensgid)
            full_transcript += get_exon_seq(exonlist[index2],ensgid)
            for i in range(index2+1,len(exonlist)):
                full_transcript += get_exon_seq(exonlist[i],ensgid)
            full_transcript_store.append(full_transcript)
        else:
            full_transcript_store.append('')
    return full_transcript_store





def detect_mode(e1,e2):
    if '_' in e1 and '_' not in e2:
        mode = 1
    elif '_' in e1 and '_' in e2:
        mode = 3
    elif '_' not in e1 and '_' in e2:
        mode = 2
    else:
        raise Exception('Can not detect mode, exons passed in second_round_match function is not Alt event type')
    return mode

def degenerate_exons(mode,e1,e2):
    if mode == 1:
        degenerated = e1.split('_')[0] + '|' + e2
    elif mode == 2:
        degenerated = e1 + '|' + e2.split('_')[0]
    elif mode == 3:
        degenerated = e1.split('_')[0] + '|' + e2.split('_')[0]
    return degenerated


def second_round_match(ensgid,exons):
    # exons is E4.1_45454-E5.6_44442
    # second round mainly solve trailing issues for Alt3 Alt5, Alt3-Alt5
    full_transcript_store = []  # ['',full_transcript1_seq,...] 
    e1,e2 = exons.split('-')
    mode = detect_mode(e1,e2)
    degenerated = degenerate_exons(mode,e1,e2)
    de1,de2 = degenerated.split('|')
    df_certain = df_exonlist.loc[df_exonlist['EnsGID']==ensgid,:]
    pattern1 = re.compile(r'\|{}\|'.format(degenerated.replace('|','\|')))  # |E3.4|E4.5|
    pattern2 = re.compile(r'^{}\|'.format(degenerated.replace('|','\|')))   # E1.1|E4.5|  first one
    pattern3 = re.compile(r'\|{}$'.format(degenerated.replace('|','\|')))  # |E6.9|E7.1  last one 
    pattern4 = re.compile(r'^{}$'.format(degenerated.replace('|','\|')))     # E4.5|E5.6   just it  
    for item in list(df_certain['Exons']): 
        if re.search(pattern1,item) or re.search(pattern2,item) or re.search(pattern3,item) or re.search(pattern4,item):
            exonlist = item.split('|')
            if mode == 1:
                transcript_left, transcript_middle, transcript_right= '', '' ,''
                for i,exon in enumerate(exonlist):
                    if exon != de1:
                        frag = get_exon_seq(exon,ensgid)
                        transcript_left += frag
                    else:
                        transcript_middle += get_trailing_exon_seq(e1,ensgid,'site1')
                        i += 1
                        for remain_exon in exonlist[i:]:
                            frag = get_exon_seq(remain_exon,ensgid)
                            transcript_right += frag
                        break
                full_transcript = transcript_left + transcript_middle + transcript_right
            elif mode == 2:
                transcript_left, transcript_middle, transcript_right= '', '' ,''
                for i,exon in enumerate(exonlist):
                    if exon != de2:
                        frag = get_exon_seq(exon,ensgid)
                        transcript_left += frag
                    else:
                        transcript_middle += get_trailing_exon_seq(e2,ensgid,'site2')
                        i += 1
                        for remain_exon in exonlist[i:]:
                            frag = get_exon_seq(remain_exon,ensgid)
                            transcript_right += frag
                        break
                full_transcript = transcript_left + transcript_middle + transcript_right    
            elif mode == 3:
                transcript_left, transcript_middle, transcript_right= '', '' ,''
                for i,exon in enumerate(exonlist):
                    if exon != de1:
                        frag = get_exon_seq(exon,ensgid)
                        transcript_left += frag
                    else:
                        transcript_middle += get_trailing_exon_seq(e1,ensgid,'site1')
                        transcript_middle += get_trailing_exon_seq(e2,ensgid,'site2')
                        i += 2
                        for remain_exon in exonlist[i:]:
                            frag = get_exon_seq(remain_exon,ensgid)
                            transcript_right += frag
                        break
                full_transcript = transcript_left + transcript_middle + transcript_right   
            full_transcript_store.append(full_transcript)
        else:
            full_transcript_store.append('')
    return full_transcript_store


def third_round_match(ensgid,exons):
    # exons should be E4.5-I4.1 or I4.5-E5.6
    # this round is for intron retention
    full_transcript_store = []  # ['',full_transcript1_seq,...] 
    e1,e2 = exons.split('-')   
    if 'I' in e1:
        i = e1; e = e2; m = 2   # m means mode, 1 means the exon is the first, 2 means the exon is the second
    elif 'I' in e2:
        i = e2; e = e1; m = 1
    df_certain = df_exonlist.loc[df_exonlist['EnsGID']==ensgid,:]
    for item in list(df_certain['Exons']):
        pattern1 = re.compile(r'{}\|'.format(e))  # E4.5|
        pattern2 = re.compile(r'^{}$'.format(e))  # E4.5$
        if re.search(pattern1,item) or re.search(pattern2,item):
            exonlist = item.split('|')  
            if m == 1:
                transcript = ''      
                for i_,exon in enumerate(exonlist):
                    if exon != e:
                        frag = get_exon_seq(exon,ensgid)
                        transcript += frag
                    else:
                        frag = get_exon_seq(e,ensgid)
                        transcript += frag
                        frag = get_exon_seq(i,ensgid)
                        transcript += frag
                        continue
            elif m == 2:
                transcript = ''      
                for i_,exon in enumerate(exonlist):
                    if exon != e:
                        frag = get_exon_seq(exon,ensgid)
                        transcript += frag
                    else:
                        frag = get_exon_seq(i,ensgid)   # fisrt add intron
                        transcript += frag
                        frag = get_exon_seq(e,ensgid)
                        transcript += frag
                        continue  
            full_transcript_store.append(transcript)      
        else:
            full_transcript_store.append('')
    return full_transcript_store


def fourth_round_match(ensgid,exons):
    # exons should be E4.5-ENSG:E4.3 or E4.5-ENSG:E4.5_4384839
    # this round is for trailing events
    full_transcript_store = []  # ['',full_transcript1_seq,...] 
    e1 = exons.split('-')[0]
    ensgid2 = exons.split('-')[1].split(':')[0]
    e2 = exons.split('-')[1].split(':')[1]
    if '_' in e2:
        e2,suffix = e2.split('_')
    else:
        suffix = None
    df_certain1 = df_exonlist.loc[df_exonlist['EnsGID']==ensgid,:]
    df_certain2 = df_exonlist.loc[df_exonlist['EnsGID']==ensgid2,:]
    for item1 in list(df_certain1['Exons']):
        for item2 in list(df_certain2['Exons']):
            pattern1_1 = re.compile(r'{}\|'.format(e1))  # E4.5|
            pattern1_2 = re.compile(r'^{}$'.format(e1))  # E4.5$
            pattern2_1 = re.compile(r'\|{}\|'.format(e2))  
            pattern2_2 = re.compile(r'^{}\|'.format(e2))  
            match_to_item1 = re.search(pattern1_1,item1) or re.search(pattern1_2,item1) 
            match_to_item2 = re.search(pattern2_1,item2) or re.search(pattern2_2,item2) 
            if match_to_item1 and match_to_item2:
                transcript = ''
                exonlist1 = item1.split('|')  
                exonlist2 = item2.split('|')
                for i1,exon1 in enumerate(exonlist1):
                    if exon1 != e1:
                        frag = get_exon_seq(exon1,ensgid)
                        transcript += frag
                    else:
                        frag = get_exon_seq(e1,ensgid)
                        for i2,exon2 in enumerate(exonlist2):
                            if exon2 == e2:
                                if suffix is not None:
                                    frag = get_trailing_exon_seq(full_e2,ensgid2,'site2')
                                else:
                                    frag = get_exon_seq(e2,ensgid2)
                                transcript += frag
                                i2 += 1
                                for remain_exon in exonlist2[i2:]:
                                    frag = get_exon_seq(remain_exon,ensgid)
                                    transcript += frag
                            else:
                                continue
                full_transcript_store.append(transcript)
            else:
                full_transcript_store.append('')
    return full_transcript_store


                



                    

