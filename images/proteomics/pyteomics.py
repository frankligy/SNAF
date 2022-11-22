#!/data/salomonis2/LabFiles/Frank-Li/proteomics_practice/pyteomics_env/bin/python3.7

from pyteomics import mzml  # mzml can be generated using ProteoWizard
from pyteomics import pylab_aux
import matplotlib.pyplot as plt
for i,dic in enumerate(mzml.read('/data/salomonis2/LabFiles/Frank-Li/neoantigen/ProteoWizard/data/OvCa114_classI_Rep#3.mzML')):

    '''
      {'index': 0, 
       'id': 'controllerType=0 controllerNumber=1 scan=1', 
       'defaultArrayLength': 6560, 
       'scanList': {'count': 1, 
                    'scan': [{'scanWindowList': {'count': 1, 
                                                 'scanWindow': [{'scan window lower limit': 400.0, 'scan window upper limit': 650.0}]}, 
                                                 'scan start time': 20.00053, 
                                                 'filter string': 'FTMS + p NSI Full ms [400.00-650.00]', 
                                                 'preset scan configuration': 1.0, 
                                                 'ion injection time': 500.0}], 'no combination': ''}, 
       'MS1 spectrum': '', 
       'ms level': 1, 
       'positive scan': '', 
       'profile spectrum': '', 
       'base peak m/z': 445.120056152344, 
       'base peak intensity': 81516.390625, 
       'total ion current': 706982.1875, 
       'lowest observed m/z': 400.00231300558, 
       'highest observed m/z': 651.566038203773,
       'count': 2, 
       'm/z array': array([400.00231301, 400.00382404, 400.00533509, ..., 651.55975553,651.56289686, 651.5660382 ]), 
       'intensity array': array([0., 0., 0., ..., 0., 0., 0.], dtype=float32)}
    '''

    if i == 8560:  # looking for msms.txt, scan index - 1
        print(dic)
        pylab_aux.plot_spectrum(dic)
        plt.savefig('test.pdf',bbox_inches='tight')
        plt.close()