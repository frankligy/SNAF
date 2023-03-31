#!/data/salomonis2/LabFiles/Frank-Li/proteomics_practice/pyteomics_env/bin/python3.7



# mzml can be generated using ProteoWizard
# ProteoWizard should be downloaded as a docker or singularity image from https://hub.docker.com/r/chambm/pwiz-skyline-i-agree-to-the-vendor-licenses
# but there are issue with singularity permission to run wine, so I either get it work using docker on my Mac or using cluster docker on /scratch
# using msconvert
# docker run -v /your/data:/data chambm/pwiz-skyline-i-agree-to-the-vendor-licenses wine msconvert file.raw

from pyteomics import mzml  
from pyteomics import pylab_aux
import matplotlib.pyplot as plt
import sys,os

original_mzml = 'VAPGEAKNL_Mel4/20140304_EXQ6_MiBa_SA_MM4-HLAp-2.mzML'
synthetic_mzml = '230323_SalomonisMix_DDA_90min_50fmol.mzML'
original_scan_id = 8443
synthetic_scan_id = 11026

for i,spectrum in enumerate(mzml.read(synthetic_mzml)):
    if i == synthetic_scan_id:
        print(spectrum)
        m_z = spectrum['m/z array'], 
        m_z = m_z[0]
        a_i = spectrum['intensity array']
        r_i = (a_i - a_i.min()) / (a_i.max() - a_i.min())
        fig,ax = plt.subplots()
        ax.bar(m_z, r_i)
        ax.set_xlabel('m/z')
        ax.set_ylabel('Relative Intensity')
        ax.set_title('MS2 Spectrum')
        plt.savefig('test.pdf',bbox_inches='tight')
        plt.close()
        break
sys.exit('stop')

# m1 = 0
# m2 = 0
# a = [(i,dic) for i,dic in enumerate(mzml.read('230323_SalomonisMix_DDA_90min_50fmol.mzML'))]
# for i,dic in a:
#     if dic['ms level'] == 1:
#         m1 += 1
#     elif dic['ms level'] == 2:
#         m2 += 1


c = 0
for i,dic in enumerate(mzml.read('230323_SalomonisMix_DDA_90min_50fmol.mzML')):
    a = dic['scanList']['scan'][0]['scan start time'].unit_info  # minute
    print(a)
    sys.exit('stop')
