__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt

ruta = '/Users/albertomolino/codigos/bpz-1.99.2/FILTER/'
ruta1 = ruta + 'GaiaDR2/'
ruta2 = ruta + 'SPLUS_July2017/'
filtros = U.get_str(ruta2+'SPLUS.list',0)
filtros_gaia = U.get_str(ruta1+'Gaiadr2.list',0)
xg1,yg1=U.get_data(ruta1+filtros_gaia[0],(0,1))
xg2,yg2=U.get_data(ruta1+filtros_gaia[1],(0,1))
xg3,yg3=U.get_data(ruta1+filtros_gaia[2],(0,1))

plt.figure(1)
plt.clf()
plt.subplot(311)
plt.plot(xg1,yg1,'b-',lw=5)
x,y = U.get_data(ruta2+filtros[0],(0,1))
plt.plot(x*1000.,y/y.max(),'-',lw=5,color='grey')
for ii in range(12):
     x,y = U.get_data(ruta2+filtros[ii],(0,1))
     plt.fill_between(x*1.,y,0,alpha=0.4,color='grey',lw=2)
     plt.plot(x*1.,y,'-',lw=1,color='grey',alpha=0.5)
     plt.xlim(3010,10900)
     #plt.xlim(3010,8300)
     plt.ylim(0.01,1.1)
     # plt.grid()
     plt.xticks(fontsize=20)
     plt.yticks(fontsize=20)
     plt.ylabel('T($\lambda$)',size=30,labelpad=5)
     # plt.xlabel('$\lambda$',size=30,labelpad=5)
     plt.legend(['Gaia-Gb','S-PLUS'],loc='upper right',fontsize=17)

plt.subplot(312)
plt.plot(xg2,yg2,'g-',lw=5)
x,y = U.get_data(ruta2+filtros[0],(0,1))
plt.plot(x*1000.,y/y.max(),'-',lw=5,color='grey')
for ii in range(12):
     x,y = U.get_data(ruta2+filtros[ii],(0,1))
     plt.plot(x*1.,y,'-',lw=1,color='grey',alpha=0.5)
     plt.fill_between(x*1.,y,0,alpha=0.4,color='grey',lw=2)
     plt.xlim(3010,10900)
     plt.ylim(0.01,1.1)
     plt.xticks(fontsize=20)
     plt.yticks(fontsize=20)
     plt.ylabel('T($\lambda$)',size=30,labelpad=5)
     # plt.xlabel('$\lambda$',size=30,labelpad=5)
     plt.legend(['Gaia-G','S-PLUS'],loc='upper right',fontsize=17)

plt.subplot(313)
plt.plot(xg3,yg3,'r-',lw=5)
x,y = U.get_data(ruta2+filtros[7],(0,1))
plt.plot(x*1000.,y/y.max(),'-',lw=5,color='grey')
for ii in range(12):
     x,y = U.get_data(ruta2+filtros[ii],(0,1))
     plt.plot(x*1.,y,'-',lw=1,color='grey',alpha=0.5)
     plt.fill_between(x*1.,y,0,alpha=0.4,color='grey',lw=2)
     #plt.xlim(5010, 10500)
     plt.xlim(3010,10900)
     plt.ylim(0.01,1.1)
     # plt.grid()
     plt.xticks(fontsize=20)
     plt.yticks(fontsize=20)
     plt.ylabel('T($\lambda$)',size=30,labelpad=5)
     plt.xlabel('$\lambda$ [$\AA$]',size=30,labelpad=5)
     plt.legend(['Gaia-Gr','S-PLUS'],loc='upper right',fontsize=17)

