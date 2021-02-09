__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import bpz_tools as B
import phz_plots as P
import matplotlib.pyplot as plt

sedlist = 'GOSMOSeB11.list'
root = '/Users/albertomolino/doctorado/articulos/SPLUS/paper0/membership/'
columns = root+'splus_auto_cali.columns'
fluxcomp = root+'members_regular.flux_comparison'
bpz = root+'members_regular.bpz'
zb,zs,chi2 = U.get_data(bpz,(1,9,8))
dz = (zb-zs)/(1.+zs)
hdf5file = root+'members_regular.hdf5'

filters = B.get_filter_list(columns)
ids = U.get_data(fluxcomp,0)
pos = N.where(ids==4539)[0][0]
z,x1,x2,x3 = P.getPDZ(hdf5file,pos)
base_y = N.arange(0.,1.,0.1)
#plt.figure(1, figsize = (6,6),dpi=75, facecolor='w', edgecolor='k')
plt.clf()
P.plotsedfittingbyID(ids[pos],columns,fluxcomp,sedlist,filters)
plt.xlabel('Wavelength $[\AA]$',size=18,labelpad=-2)
plt.ylabel('Magnitude',size=20,labelpad=2)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
#plt.ylim(23.,17)
plt.xlim(3200,9300)
siete = plt.axes([0.5,.215,0.375,0.375])
y = x3/x3.sum()
plt.plot(z,y/y.sum(),'-',color='black',linewidth=1,alpha=0.99)
plt.fill_between(z,y/y.sum(),0.,color='grey',alpha=0.5)
plt.ylabel('p(z)',size=20,labelpad=15)
plt.xlabel('$z$',size=30,labelpad=0)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylim(0.0001,0.0249)
plt.grid()
plt.xlim(0.01,0.199)
plt.plot(base_y*0.+0.050,base_y,'r--',lw=2)
label_1 = '$z_{b}=%.3f$'%(zb[pos])
label_2 = '$z_{s}=%.3f$'%(zs[pos])
label_3 = '$dz/1+z=%.1f$'%(dz[pos]*100.)+'%'
label_4 = '$\chi^{2}=%.3f$'%(chi2[pos])
plt.legend([label_1+'\n'+label_2+'\n'+label_3+'\n'+label_4,'$z_{cluster}=0.050$'],loc='upper right',fontsize=11)
#plt.title('clash_a209_nir_0452',size=25)
plt.savefig('/Users/albertomolino/Desktop/SED%i.png'%(ids[pos]),dpi=100)