__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
#import bpz_tools as B

# Paths
root = '/Volumes/CLASH/S82/specz/AbsMag/analysis/'
catalog = root + 'specz_vs_photoz_MassMr.cat'

# Colors for data
scolors = list(cm.RdBu_r(N.linspace(0, 1, 12)))

# Reading data
Mr_sp,Mr_ph,Mas_sp,Mas_ph,zs,SpT = U.get_data(catalog,(2,8,3,9,4,5))


## Plots
# Stellar_Mass
plt.figure(1,figsize = (8.5,6.7),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
plt.scatter(zs[::3]-100.,Mas_sp[::3],s=100,c=SpT[::3],marker=u'.',cmap=cm.RdBu,alpha=0.5)
cb = plt.colorbar(pad=0.,format='%i',ticks=N.arange(1.,17,1))
cb.set_label(label='Spectral-type',size=20)
#cb.yticks(fontsize=20)
plt.scatter(zs,Mas_sp,s=150,c=SpT,marker=u'+',cmap=cm.RdBu,alpha=0.25)
plt.xlim(0.,0.8)
plt.ylim(8.01,11.9)
plt.grid()
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('$z$',size=30,labelpad=-5)
plt.ylabel('Stellar Mass log M$_{*}$ ($M_{\odot}$)',size=20,labelpad=5)
#cb.set_label('Integrated-PDF [$\%$]',size=28,labelpad=20)

plt.figure(2,figsize = (8.5,6.7),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
basem = N.arange(-2,2,0.01)
dmass = (Mas_sp-Mas_ph)
plt.clf()
a1,a2,a3 = plt.hist(dmass+0.041,basem,facecolor='blue',alpha=0.25,normed=1)
plt.yticks([])
plt.annotate('M$_{*,specz}$ - M$_{*,photoz}$',xy=(-1.45,3.1),fontsize=40,color='purple')
plt.xticks(N.arange(-2,3,1),('-2.0','-1.0','0.0','1.0','2.0'),size=24)


## Plots
# Absolute_Mags
plt.figure(3,figsize = (8.5,6.7),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
plt.scatter(zs[::3]-100.,Mr_sp[::3],s=100,c=SpT[::3],marker=u'.',cmap=cm.RdBu,alpha=0.5)
cb = plt.colorbar(pad=0.,format='%i',ticks=N.arange(1.,17,1))
cb.set_label(label='Spectral-type',size=20)
#cb.yticks(fontsize=20)
plt.scatter(zs,Mr_sp,s=150,c=SpT,marker=u'+',cmap=cm.RdBu,alpha=0.25)
plt.xlim(0.,0.8)
plt.ylim(-16.1,-24.9)
plt.grid()
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('$z$',size=30,labelpad=-5)
plt.ylabel('$M_{r}$',size=24,labelpad=5)
#cb.set_label('Integrated-PDF [$\%$]',size=28,labelpad=20)

plt.figure(4,figsize = (8.5,6.7),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
basem = N.arange(-3,3,0.03)
dmr = (Mr_sp-Mr_ph)
a1,a2,a3 = plt.hist(dmr+0.0,basem,facecolor='blue',alpha=0.25,normed=1)
plt.yticks([])
plt.ylim(0.,1.3)
plt.annotate('M$_{r,specz}$ - M$_{r,photoz}$',xy=(-2.1,1.15),fontsize=40,color='purple')
plt.xticks(N.arange(-2,3,1),('-2.0','-1.0','0.0','1.0','2.0'),size=24)
plt.savefig('/Users/albertomolino/Desktop/Mrzs2.png',dpi=100)



