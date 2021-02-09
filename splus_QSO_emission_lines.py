__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

root_to_bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root_to_filts = root_to_bpz+'FILTER/'
filters = U.get_str(root_to_filts+'splusNBs.list',0)

narrow = []
narrow.append('F378')
narrow.append('F395')
narrow.append('F410')
narrow.append('F430')
narrow.append('F515')
narrow.append('F660')
narrow.append('F861')

#QSOs
galq = N.zeros(10)
galq[0]=1034
galq[1]=1215
galq[2]=1240
galq[3]=1549
galq[4]=1908
galq[5]=2799
galq[6]=4341
galq[7]=4862
galq[8]=5008
galq[9]=6564

lab = ['OVI: 1034 $\AA$','Ly-alpha: 1215 $\AA$','NV: 1240 $\AA$',
 'CIV: 1549 $\AA$','CIII: 1908 $\AA$','MgII: 2799 $\AA$','H-gama: 4341 $\AA$',
 'H-beta: 4863 $\AA$','OIII: 5008 $\AA$','H-alpha: 6564 $\AA$']

# Redshift
z = N.arange(0.,10.01,0.01)

colores = list(cm.jet(N.linspace(0, 1, 10)))
colores = colores[::-1]

plt.figure(1, figsize=(11,9),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
for ii in range(10):
    plt.plot(galq[ii]*(1+z),z,'-',alpha=0.5,lw=3.,color=colores[ii])

for ii in range(7):
        x,y = U.get_data(root_to_filts+filters[ii],(0,1))
        plt.plot(x,8.*(y/y.max()),'-',lw=1,color='grey',alpha=0.5)
        plt.fill_betweenx(8.*(y/y.max()),x,0,alpha=0.1)
        plt.xlabel('Wavelength [$\AA$]',size=25,labelpad=5)
        plt.ylabel('$z$',size=35)

for ii in range(10):
    plt.plot(galq[ii]*(1+z),z,'-',alpha=0.5,lw=3,color=colores[ii])
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(3500,13000)
plt.ylim(1.,10.)
plt.legend(lab,loc='upper right',fontsize=30)
plt.grid()
