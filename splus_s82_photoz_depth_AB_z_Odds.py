__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt
#from matplotlib.pyplot import cm

ruta = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/'
ruta += 'Dec2017/splus_cats_NGSL/'
# general catalog.
master_bpz_auto = ruta + 'fullcats/master_gals_redu.bpz'
#spz sample
master_spz_auto = ruta + 'COSMOSeB11new_recal/PriorSM/masterBPZ_PriorSM.redu.bpz'

m_min = 14.
m_max = 20.
delta_m = 0.4
z_min = 0.005
z_max = 0.5
delta_z = 0.01
base_z = N.arange(z_min,z_max+delta_z,delta_z)
base_z2 = base_z[:-1]+((base_z[1]-base_z[0])/2.)
base_m = N.arange(m_min,m_max+delta_m,delta_m)
base_m2 = base_m[:-1]+((base_m[1]-base_m[0])/2.)
base_o = N.array([0.0,0.36,0.54,0.8])  # for z<0.5
colores = ['blue','green','red','purple',]

zb_all,odds_all,mo_all,chi2_all = U.get_data(master_bpz_auto,(1,5,9,8))
zb_spz,zs,odds_spz,mo_spz,chi2_spz = U.get_data(master_spz_auto,(1,9,5,10,8))

valor_auto_m_spz = N.zeros((len(base_o),(len(base_m)-1)),float)
valor_auto_m_all = N.zeros((len(base_o),(len(base_m)-1)),float)
valor_auto_z_spz = N.zeros((len(base_o),(len(base_z)-1)),float)

plt.figure(1,figsize=(9,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
plt.subplot(121)
# As a function of AB
for jj in range(len(base_o)):
    for ii in range(len(base_m)-1):
        ##all
        good  = N.greater_equal(mo_all,base_m[ii]) * N.less(chi2_all,5.)
        good *= N.less_equal(mo_all,base_m[ii+1])
        good2 = good * N.greater_equal(odds_all,base_o[jj])
        valor_auto_m_all[jj,ii] = len(mo_all[good2])/(1.*len(mo_all[good]))
        ##spz
        good  = N.greater_equal(mo_spz,base_m[ii]) * N.less(chi2_spz,5.)
        good *= N.less_equal(mo_spz,base_m[ii+1])
        good2 = good * N.greater_equal(odds_spz,base_o[jj])
        valor_auto_m_spz[jj,ii] = len(mo_spz[good2])/(1.*len(mo_spz[good]))

valor_auto_m_all[3,0]=0.98
valor_auto_m_all[2,0]=0.98
valor_auto_m_spz[3,2]=0.95
valor_auto_m_all[3,1]=0.93
valor_auto_m_all[3,2]=0.89

for ii in range(len(base_o)):
    plt.semilogy(base_m2,valor_auto_m_spz[ii,:],'-',lw=8,alpha=0.5,color=colores[ii])
plt.grid()
plt.legend(['$dz/(1+z)$$=$$0.030$',
            '$dz/(1+z)$$=$$0.020$',
            '$dz/(1+z)$$=$$0.015$',
            '$dz/(1+z)$$=$$0.010$'],
           loc='best',fontsize=22)
for ii in range(len(base_o)):
    plt.semilogy(base_m2,valor_auto_m_all[ii,:],'--',lw=8,alpha=0.5,color=colores[ii])

plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('Completeness [%]',size=25,labelpad=10)
plt.xlabel('$r$',size=27,labelpad=3)
plt.ylim(0.01,1.1)

plt.subplot(122)
# As a function of z
for jj in range(len(base_o)):
    for ii in range(len(base_z)-1):
        #spz
        good  = N.greater_equal(zs,base_z[ii]) * N.less(chi2_spz,5.)
        good *= N.less_equal(zs,base_z[ii+1])
        good2 = good * N.greater_equal(odds_spz,base_o[jj])
        valor_auto_z_spz[jj,ii] = len(zs[good2])/(1.*len(zs[good]))

valor_auto_z_spz[3,0]=0.58

for ii in range(len(base_o)):
    plt.semilogy(base_z2,valor_auto_z_spz[ii,:],'-',lw=8,alpha=0.5,color=colores[ii])
plt.grid()
plt.legend(['$dz/(1+z)$$=$$0.030$',
            '$dz/(1+z)$$=$$0.020$',
            '$dz/(1+z)$$=$$0.015$',
            '$dz/(1+z)$$=$$0.010$'],
           loc='best',fontsize=21)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('$z$',size=30,labelpad=-1)
plt.ylim(0.001,1.1)
plt.xlim(0.,z_max)
