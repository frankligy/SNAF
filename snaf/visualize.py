#!/data/salomonis2/LabFiles/Frank-Li/refactor/neo_env/bin/python3.7

import numpy as np
import pandas as pd
import os
import sys
import pickle
import h5py
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle 
import matplotlib as mpl

'''
contain visualization functions
'''

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'

def show_candicates(ax,aa,extra,n_from_first,hla,phase,evidences,first,second,dna_first,dna_second,binding_score,immunogenicity_score):
    # set lim
    ax.set_xlim(0,100)
    ax.set_ylim(0,100)
    # draw exons
    rect_first = Rectangle((30,10),20,5,linewidth=1,edgecolor='k',facecolor='blue')
    rect_second = Rectangle((50,10),20,5,linewidth=1,edgecolor='k',facecolor='orange')
    ax.add_patch(rect_first)
    ax.add_patch(rect_second)
    # draw junction seq
    seq_first_text = ax.text(x=50,y=20,s=first[-dna_first:],color='blue',va='bottom',ha='right',fontsize=2)
    seq_second_text = ax.text(x=50,y=20,s=second[:dna_second],color='orange',va='bottom',ha='left',fontsize=2)
    plt.pause(0.01)
    # draw aa seq
    bbox_coords = ax.transData.inverted().transform(seq_first_text.get_window_extent().get_points())
    width = bbox_coords[1,0] - bbox_coords[0,0]
    height = bbox_coords[1,1] - bbox_coords[0,1]
    start_x = bbox_coords[0,0] + width / (2 * dna_first) * 3
    start_y = 20 + height + 5
    tmp_aa_list = []
    tmp_aa_list[:] = aa
    aa_to_draw = '  '.join(tmp_aa_list)
    aa_text = ax.text(x=start_x,y=start_y,s=aa_to_draw,color='r',fontweight='bold',fontsize=2)
    plt.pause(0.01)
    # draw barplot for scores
    barcontainer = ax.bar(x=[25,75],height=[np.interp(binding_score,xp=(0,2),fp=(0,50)),np.interp(immunogenicity_score,xp=(0,1),fp=(0,50))],bottom=50,width=20,color=['#158BFB','#F2075D'])
    ax.text(x=barcontainer[0].get_x() + barcontainer[0].get_width()/2,y=50,s='binding',fontsize=3,va='top',ha='center')
    ax.text(x=barcontainer[1].get_x() + barcontainer[1].get_width()/2,y=50,s='immunogenicity',fontsize=3,va='top',ha='center')
    ax.text(x=barcontainer[0].get_x() + barcontainer[0].get_width()/2,y=barcontainer[0].get_y() + barcontainer[0].get_height(),s='{}'.format(round(binding_score,3)),fontsize=3,va='top',ha='center')
    ax.text(x=barcontainer[1].get_x() + barcontainer[1].get_width()/2,y=barcontainer[1].get_y() + barcontainer[1].get_height(),s='{}'.format(round(immunogenicity_score,3)),fontsize=3,va='top',ha='center')
    # annotate HLA and score
    ax.set_title('{}\n{}\n{}'.format(hla,phase,evidences),fontsize=3)
    # remove tick and labels
    ax.set_xticks([])
    ax.set_yticks([])
    return ax

def get_base_subexon_and_trail(subexon):
    if 'U' in subexon:
        post = subexon
        trail = None
    elif 'ENSG' in subexon:
        post = subexon
        trail = None
    elif '_' in subexon:
        post,trail = subexon.split('_')
    else:
        post = subexon
        trail = None
    return post,trail
    


def draw_genome(ax,uid,dict_exonCoords):
    ensid = uid.split(':')[0]
    first,second = uid.split(':')[1].split('-')
    first,trail1 = get_base_subexon_and_trail(first)
    second,trail2 = get_base_subexon_and_trail(second)
    df_all_subexons = pd.DataFrame(data=dict_exonCoords[ensid]).T
    df_all_subexons.columns = ['chr','strand','start','end','suffer']
    df_all_subexons.sort_values(by='start',inplace=True)
    strand = df_all_subexons.iloc[0,1]
    chr_ = df_all_subexons.iloc[0,0]
    # adjust the ax
    ax.set_xlim(-0.05,1.05)
    ax.set_ylim(-0.05,1.05)
    if strand == '+':
        starting, ending = df_all_subexons.iloc[0,2], df_all_subexons.iloc[-1,3]
        for row in df_all_subexons.itertuples():
            subexon_start = np.interp(x=row.start,xp=[starting,ending],fp=[0,1])
            subexon_end = np.interp(x=row.end,xp=[starting,ending],fp=[0,1])
            if row.Index == first or row.Index == second:
                facecolor = 'red'
                edgecolor = 'k'
            else:
                facecolor = '#0D1273'
                edgecolor = 'k'
            if row.Index.startswith('E'):
                subexon_rect = Rectangle((subexon_start,0.3),subexon_end - subexon_start, 0.4,linewidth=0.1,facecolor=facecolor,edgecolor=edgecolor)
                ax.add_patch(subexon_rect)
                ax.text(subexon_start,0.3,row.Index,fontsize=1,va='top',ha='left',rotation=60) 
            elif row.Index.startswith('I'):
                subexon_rect = Rectangle((subexon_start,0.45),subexon_end - subexon_start, 0.1,linewidth=0.1,facecolor=facecolor,edgecolor=edgecolor)
                ax.add_patch(subexon_rect)
                ax.text(subexon_start,0.3 + 0.4,row.Index,fontsize=1,va='bottom',ha='left',rotation=60)        

    else:
        starting, ending = df_all_subexons.iloc[0,3], df_all_subexons.iloc[-1,2] 
        for row in df_all_subexons.itertuples():
            subexon_start = np.interp(x=row.end,xp=[starting,ending],fp=[0,1])
            subexon_end = np.interp(x=row.start,xp=[starting,ending],fp=[0,1])
            if row.Index == first or row.Index == second:
                facecolor = 'red'
                edgecolor = 'k'
            else:
                facecolor = '#0D1273'
                edgecolor = 'k'
            if row.Index.startswith('E'):
                subexon_rect = Rectangle((subexon_start,0.3),subexon_end - subexon_start, 0.4,linewidth=0.1,facecolor=facecolor,edgecolor=edgecolor)
                ax.add_patch(subexon_rect)
                ax.text(subexon_start,0.3,row.Index,fontsize=1,va='top',ha='right',rotation=60) 
            elif row.Index.startswith('I'):
                subexon_rect = Rectangle((subexon_start,0.45),subexon_end - subexon_start, 0.1,linewidth=0.1,facecolor=facecolor,edgecolor=edgecolor)
                ax.add_patch(subexon_rect)
                ax.text(subexon_start,0.3 + 0.4,row.Index,fontsize=1,va='bottom',ha='right',rotation=60)     

    # final touch
    ax.set_xticks([0,1])
    ax.set_xticklabels([starting,ending])     
    ax.text(0.5,1,chr_,va='top',ha='center',fontsize=5)  
    ax.set_yticks([])
    for trail in [trail1,trail2]:
        if trail is not None:
            ax.axvline(x=np.interp(x=int(trail),xp=[starting,ending],fp=[0,1]),ymin=0.3,ymax=0.7,linewidth=0.1,linestyle='--')
            ax.text(x=np.interp(x=int(trail),xp=[starting,ending],fp=[0,1]),y=0.2,s=trail,va='top',ha='center',fontsize=1,rotation=60)            
    return ax




















