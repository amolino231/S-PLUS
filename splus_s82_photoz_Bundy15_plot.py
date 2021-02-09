__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt

cat = '/Volumes/FLASHDRIVE/matched_spzsample.cat'
mr,zs,zb,zb82 = U.get_data(cat,(4,6,5,7))
dz_82 = (zb82-zs)/(1.+zs)
dz_s = (zb-zs)/(1.+zs)

res=0.004
basez = N.arange(-2.,2.,res)
basez2 = basez[:-1]+((basez[1]-basez[0])/2.)
b1,b2,b3 = plt.hist(dz_82[mr<19],basez,color='blue',alpha=0.25)
a1,a2,a3 = plt.hist(dz_s[mr<19],basez,color='red',alpha=0.25)
plt.plot(basez2,b1,'b-',basez2,a1,'r-',lw=3)
plt.xlim(-0.15,0.15)
plt.grid()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel('$\delta_{z}/(1+z)$',size=28,labelpad=3)
plt.ylabel('$\#$',size=28)
label1=' ''\n''$\sigma_{z}=0.031$''\n''$\mu_{z}=0.027$'
label2=' ''\n''$\sigma_{z}=0.016$''\n''$\mu_{z}=0.000$'
plt.legend([label1,label2],loc='upper left',fontsize=25)