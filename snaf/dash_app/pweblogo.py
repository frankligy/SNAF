import matplotlib as mpl
mpl.use('Agg')
from matplotlib.patches import PathPatch
import matplotlib.pyplot as plt
from scipy.stats import entropy
from collections import Counter
import seaborn as sns
import numpy as np
import argparse
from .schema import *
import os,sys

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


def draw_letter(letter,x,y,x_scale=1,y_scale=1,ax=None):
    text = letter_schema[letter]
    transform = mpl.transforms.Affine2D().scale(x_scale*globscale,y_scale*globscale) + \
                mpl.transforms.Affine2D().translate(x,y) + \
                ax.transData
    p = PathPatch(text,lw=0,fc=color_schema[letter],transform=transform)
    ax.add_artist(p)

def draw_weblogo(all_score):
    fig,ax = plt.subplots()
    x = 1
    highest = 0
    for score in all_score:
        y = 0
        for b,s in score:
            draw_letter(b,x,y,x_scale=1,y_scale=s,ax=ax)
            #raise Exception
            y += s
        x += 1.5
        highest = max(highest,y)
    ax.set_xlim((0,x))
    ax.set_ylim((0,highest))
    ax.set_ylabel('bits')
    ax.set_xticks(np.linspace(1,1+(len(all_score)-1)*1.5,len(all_score)))
    ax.set_xticklabels(['pos{}'.format(i) for i in range(1,len(all_score)+1)])
    ax.set_xlabel('position')
    sns.despine(ax=ax,offset=10,trim=True)

def calculate_weblogo(pep_list):
    amino = 'ARNDCQEGHILKMFPSTWYV'
    n = len(pep_list)
    l = len(pep_list[0])
    all_score = []
    all_bit = []
    for i,pos in enumerate(zip(*pep_list)):
        score = []
        counter = Counter(pos)
        freq = np.array(list(counter.values()))/n
        ent = np.log2(20) - entropy(freq,base=2)
        for aa in amino:
            score.append((aa,ent*counter.get(aa,0)/n))
        score = sorted(score,key=lambda x:x[1])
        all_score.append(score)
        all_bit.append(ent)
    return all_score,all_bit

def run_pweblogo(pep_list,out):
    all_score,all_bit = calculate_weblogo(pep_list)
    draw_weblogo(all_score)
    plt.savefig(out,bbox_inches='tight')
    plt.close()


