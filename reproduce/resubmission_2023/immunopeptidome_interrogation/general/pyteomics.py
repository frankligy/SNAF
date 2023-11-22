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

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'


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
    elif ion_type == 'M':
        a_seq = seq
        m_z = mass.calculate_mass(sequence=a_seq,ion_type=ion_type,charge=charge)
        color = 'orange'
        assemble = ion_type + '+' * charge
    return m_z,color, assemble

def cosine_similarity(a,b):
    dot = np.dot(a, b)
    norma = np.linalg.norm(a)
    normb = np.linalg.norm(b)
    cos = dot / (norma * normb)
    return cos


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
        axes[1].text(x=m_z[idx],y=r_i[idx]-0.05,s=assemble,va='top',ha='center',fontsize=4,c=color)

    # calculate smilarity
    p, s = pearsonr(x=store_intensity_exp,y=store_intensity_syn)
    cos = cosine_similarity(store_intensity_exp,store_intensity_syn)
    fig.suptitle('{}\nPearson r:{}\np-value:{}\ncosine:{}'.format(seq,round(p,2),round(s,2),cos),fontsize=6)
    
    plt.savefig('mirror_{}.pdf'.format(seq),bbox_inches='tight')
    plt.close()
            




# VAPGEAKNL_Mel4
exp_mzml = 'VAPGEAKNL_Mel4/20140304_EXQ6_MiBa_SA_MM4-HLAp-2.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 8443
syn_scan_num = 11026
mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,800),[('y',1,1),('y',2,1),('y',3,1),('y',6,1),('y',7,1),('y',7,2),('y',8,2),('b',2,1)],'VAPGEAKNL')


# YALANIKWI_Mel12
exp_mzml = 'YALANIKWI_Mel12/20140304_EXQ6_MiBa_SA_MM12-HLAp-2.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 15779
syn_scan_num = 31110
mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,1000),[('y',1,1),('y',6,1),('y',7,1),('y',8,1),('b',2,1),('b',3,1),('b',8,1),('M',None,2)],'YALANIKWI')


# HAAASFETL_Mel15, common
exp_mzml = 'HAAASFETL_Mel15/20141208_QEp7_MiBa_SA_HLA-I-p_MM15_3_A.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 35699
syn_scan_num = 18476
mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,900),[('y',1,1),('y',2,1),('y',3,1),('b',2,1),('b',3,1),('b',4,1),('b',5,1),('b',6,1),('b',7,1),('b',8,1)],'HAAASFETL')

# TELQRTLSL_Mel26
exp_mzml = 'TELQRTLSL_Mel26/20141216_QEp7_MiBa_SA_HLA-I-p_MMf_3_1.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 33048
syn_scan_num = 20400
mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,1000),[('y',1,1),('y',2,1),('y',5,1),('y',6,1),('y',7,1),('y',8,1),('b',2,1),('b',3,1),('b',5,1),('b',6,1),('b',7,1),('b',8,1)],'TELQRTLSL')

# IIVDQKQLV_Mel25
exp_mzml = 'IIVDQKQLV_Mel25/20141216_QEp7_MiBa_SA_HLA-I-p_MMf_2_1.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 31896
syn_scan_num = 19422
mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,1000),[('y',1,1),('y',2,1),('y',7,1),('y',8,1),('b',2,1),('b',6,1),('b',7,1),('b',8,1)],'IIVDQKQLV')

# AEVPYRVDL_Mel26
exp_mzml = 'AEVPYRVDL_Mel26/20141216_QEp7_MiBa_SA_HLA-I-p_MMf_3_2.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 36153
syn_scan_num = 21468
mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,1000),[('y',1,1),('y',6,1),('y',7,1),('b',2,1),('b',3,1),('b',8,1)],'AEVPYRVDL')

# KEKLDQLVY_Mel27
exp_mzml = 'KEKLDQLVY_Mel27/20141214_QEp7_MiBa_SA_HLA-I-p_MMf_4_2.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 22277
syn_scan_num = 18444
mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,1000),[('y',1,1),('y',2,1),('b',2,1),('b',3,1),('b',4,1),('b',6,1),('b',7,1),('b',8,1)],'KEKLDQLVY')

# TLPVQPAEL_Mel39
exp_mzml = 'TLPVQPAEL_Mel39/20141215_QEp7_MiBa_SA_HLA-I-p_MMf_16_1.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 47073
syn_scan_num = 22882
mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,900),[('y',1,1),('y',2,1),('y',4,1),('y',5,1),('y',7,1),('b',2,1),('b',3,1),('b',4,1),('b',6,1),('b',7,1),('b',8,1)],'TLPVQPAEL')

# VVPHACNASY_OvCa111
exp_mzml = 'VVPHACNASY_OvCa111/OvCa111_classI_Rep#4.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 11630
syn_scan_num = 12961
mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,1000),[('y',1,1),('y',2,1),('y',3,1),('y',8,1),('y',9,1),('b',2,1),('b',6,1),('b',8,1),('b',9,1)],'VVPHACNASY')

# KGPWYPLSL_OvCa80
exp_mzml = 'KGPWYPLSL_OvCa80/OvCa80_classI_Rep#2.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 13509
syn_scan_num = 30214
mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,1000),[('y',2,1),('y',4,1),('b',5,1),('b',6,1),('b',7,1),('b',8,1),('b',8,2),('b',6,2),('b',5,2)],'KGPWYPLSL')

# NQDEDPLEV_OvCa99
exp_mzml = 'NQDEDPLEV_OvCa99/OvCa99_classI_Rep#1.mzML'
syn_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
exp_scan_num = 13051
syn_scan_num = 20662
mirror_plot(exp_mzml,syn_mzml,exp_scan_num,syn_scan_num,(0,1000),[('y',2,1),('y',3,1),('y',4,1),('y',5,1),('b',2,1),('b',6,1),('b',7,1),('b',8,1)],'NQDEDPLEV')

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

