import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('/Volumes/GoogleDrive/My Drive/SNAF/revision_STM/SNAF_B/SM FL Neo/MET Ref 1/Book1.csv',index_col=0)
# plt.style.use('dark_background')

fig,ax = plt.subplots()
d = data['Distance [µm]'].tolist()
n = data['DAPI 512(1)'].tolist()
s = data['GFP 512(1)']
c = data['Cy5 512(1)']
ax.plot(d,n,color='#0C035C',marker='',linestyle='-')
ax.plot(d,s,c='#3D813E',marker='',linestyle='-')
ax.plot(d,c,c='#74101D',marker='',linestyle='-')
ax.set_ylim([0,12000])



