#!/data/salomonis2/LabFiles/Frank-Li/proteomics_practice/pyteomics_env/bin/python3.7

# mzml can be generated using ProteoWizard
# ProteoWizard should be downloaded as a docker or singularity image from https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
# but there are issue with singularity permission to run wine, so I either get it work using docker on my Mac or using cluster docker on /scratch
# using msconvert
# docker run -v /your/data:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert file.raw

from pyteomics import mzml  
from pyteomics import mass
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr
import sys,os


# functions
def calculate_m_z(ion,seq):
    ion_type, position, charge = ion
    if ion_type == 'y':
        a_seq = seq[-position:]
        m_z = mass.calculate_mass(sequence=a_seq,ion_type=ion_type,charge=charge)
        color = 'r'
        assemble = ion_type + str(position) + '+' * charge
    elif ion_type == 'b':
        a_seq = seq[:position]
        m_z = mass.calculate_mass(sequence=a_seq,ion_type=ion_type,charge=charge)
        color = 'b'
        assemble = ion_type + str(position) + '+' * charge
    return m_z,color, assemble




def mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,xlim,ions,seq):
    fig,axes = plt.subplots(nrows=2,ncols=1,sharex=True,gridspec_kw={'hspace':0.1,'height_ratios':(0.5,0.5)})
    store_intensity_exp = []
    store_intensity_syn = []
    # first draw exp
    for i,spectrum in enumerate(mzml.read(exp_mzml)):
        if i == exp_scan_num:
            m_z = spectrum['m/z array']
            a_i = spectrum['intensity array']
            r_i = (a_i - a_i.min()) / (a_i.max() - a_i.min())
            axes[0].bar(m_z,r_i,color='k')
            axes[0].set_ylabel('Relative Intensity Exp')
            axes[0].spines['top'].set_visible(False)
            axes[0].spines['right'].set_visible(False)
            axes[0].tick_params(axis='x',length=0)
            axes[0].set_xlim(xlim)
            break
        
    for ion in ions:
        theo_m_z, color, assemble = calculate_m_z(ion,seq)
        idx = np.abs((m_z - theo_m_z)).argmin()
        store_intensity_exp.append(r_i[idx])
        axes[0].vlines(x=m_z[idx],ymin=0,ymax=r_i[idx],colors=[color],linewidth=1.5)
        axes[0].text(x=m_z[idx],y=r_i[idx]+0.05,s=assemble,va='bottom',ha='center',fontsize=4,c=color)


    # then draw syn
    for i,spectrum in enumerate(mzml.read(syn_mzml)):
        if i == syn_scan_num:
            m_z = spectrum['m/z array']
            a_i = spectrum['intensity array']
            r_i = np.negative((a_i - a_i.min()) / (a_i.max() - a_i.min()))
            axes[1].spines['top'].set_position(('data',0))
            axes[1].bar(m_z,r_i,color='k')
            axes[1].set_ylabel('Relative Intensity Syn')
            axes[1].spines['bottom'].set_visible(False)
            axes[1].spines['right'].set_visible(False)
            break

    for ion in ions:
        theo_m_z, color, assemble = calculate_m_z(ion,seq)
        idx = np.abs((m_z - theo_m_z)).argmin()
        store_intensity_syn.append(-r_i[idx])
        axes[1].vlines(x=m_z[idx],ymin=0,ymax=r_i[idx],colors=[color],linewidth=1.5)
        axes[1].text(x=m_z[idx],y=r_i[idx]-0.05,s=assemble,va='top',ha='center',fontsize=8,c=color)

    # calculate pearsonr
    p, s = pearsonr(x=store_intensity_exp,y=store_intensity_syn)
    fig.suptitle('{}\nPearson r:{}\np-value:{}'.format(seq,round(p,2),round(s,2)),fontsize=4)

    plt.savefig('test.pdf',bbox_inches='tight')
    plt.close()
            


exp_mzml = 'VAPGEAKNL_Mel4/20140304_EXQ6_MiBa_SA_MM4-HLAp-2.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 8443
syn_scan_num = 11026

mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,800),[('y',1,1),('y',2,1),('y',3,1),('y',6,1),('y',7,1),('y',7,2),('y',8,2),('b',2,1)],'VAPGEAKNL')
sys.exit('stop')






# calculate fragment ion m/z
from pyteomics import mass
sequence = 'VAPGEAKNL'
sequence = 'PGEAKNL'
ion_type = 'y'
charge = 1
a = mass.Composition(sequence='APGEAKNL')
b = mass.Composition(formula='H2O')
m_z = mass.calculate_mass(composition=a-b,ion_type='y',charge=1)
m_z = mass.calculate_mass(sequence='VA',ion_type='a',charge=charge)


for i,dic in enumerate(mzml.read('230323_SalomonisMix_DDA_90min_50fmol.mzML')):
    a = dic['scanList']['scan'][0]['scan start time'].unit_info  # minute

