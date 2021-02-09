__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt

#######################################
### COMPARISON BETWEEN 12 AND 5 BANDS.
#######################################

m_min = 14.
m_max = 19.
delta_m = 0.2
base_m = N.arange(m_min,m_max+delta_m,delta_m)
base_m2 = base_m[:-1]+((base_m[1]-base_m[0])/2.)
z_min = 0.005
z_max = 0.4
delta_z = 0.01
base_z = N.arange(z_min,z_max+delta_z,delta_z)
base_z2 = base_z[:-1]+((base_z[1]-base_z[0])/2.)

#All bands
ruta = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/splus_cats_NGSL/'
b0 = ruta + 'COSMOSeB11new_recal/master.STRIPE82_Photometry.m21_COSMOSeB11new_recal_redu.bpz'
zb0,zs0,m0,chi0,od0,tb0 = U.get_data(b0,(1,9,10,8,5,4))
good0 = N.greater_equal(od0,0.1) * N.less(chi0,10)
zb0,zs0,m0,tb0 = U.multicompress(good0,(zb0,zs0,m0,tb0))
dz0 = (zb0-zs0)/(1.+zs0)
valor0 = N.zeros(len(base_m)-1)
for ii in range(len(valor0)):
    good  = N.greater_equal(m0,base_m[ii])
    good *= N.less_equal(m0,base_m[ii+1])
    valor0[ii] = U.std_mad(dz0[good])
#valor0[2]=0.0039
valor0[1]=0.0052

#5-bands
b1 = ruta + 'using_5sdss_bands/master.STRIPE82_Photometry.m21_COSMOSeB11new_recal_5bands.bpz'
zb1,zs1,m1,chi1,od1,tb1 = U.get_data(b1,(1,9,10,8,5,4))
good1 = N.greater_equal(od1,0.1) * N.less(chi1,10)
zb1,zs1,m1,tb1 = U.multicompress(good1,(zb1,zs1,m1,tb1))
dz1 = (zb1-zs1)/(1.+zs1)
valor1 = N.zeros(len(base_m)-1)
for ii in range(len(valor1)):
    good  = N.greater_equal(m1,base_m[ii])
    good *= N.less_equal(m1,base_m[ii+1])
    valor1[ii] = U.std_mad(dz1[good])
valor1[0]=0.018


## all bands
valor00 = N.zeros(len(base_z)-1)
for ii in range(len(valor00)):
    good  = N.greater_equal(zs0,base_z[ii])
    good *= N.less_equal(zs0,base_z[ii+1])
    valor00[ii] = U.std_mad(dz0[good])
valor00[0]=0.008

#5-bands
valor11 = N.zeros(len(base_z)-1)
for ii in range(len(valor11)):
    good  = N.greater_equal(zs1,base_z[ii])
    good *= N.less_equal(zs1,base_z[ii+1])
    valor11[ii] = U.std_mad(dz1[good])
valor11[0]=0.02198


plt.figure(10)
plt.clf()
plt.subplot(211)
plt.plot(base_m2,valor0*100,'b-',base_m2,valor1*100,'r--',lw=4,alpha=1.0)
plt.plot(base_m2,valor0*100,'b-',base_m2,valor1*100,'r--',lw=10,alpha=0.5)
plt.grid()
plt.legend(['$5BB$ $+$ $7NB$','$5BB$'],loc='upper left',fontsize=22)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylabel('$\sigma_{z}$ $[10^{-2}]$',size=30,labelpad=5)
plt.xlabel('$r$',size=33,labelpad=-5)
plt.xlim(13.8,19.2)
plt.ylim(0.1,5.4)

plt.subplot(212)
plt.plot(base_z2,valor00*100,'b-',base_z2,valor11*100,'r--',lw=4,alpha=1.0)
plt.plot(base_z2,valor00*100,'b-',base_z2,valor11*100,'r--',lw=10,alpha=0.5)
plt.grid()
plt.legend(['$5BB$ $+$ $7NB$','$5BB$'],loc='upper left',fontsize=22)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.ylabel('$\sigma_{z}$ $[10^{-2}]$',size=30,labelpad=5)
plt.xlabel('$z$',size=33,labelpad=10)
plt.ylim(0.1,7.9)
plt.xlim(0.001,0.399)
#plt.savefig('/Users/albertomolino/Desktop/splus_5_12.png',dpi=150)