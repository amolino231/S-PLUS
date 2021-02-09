__author__ = 'albertomolino'

## This code is stored at: ~/Postdoc/T80S_Pipeline/Commisioning/codes/

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt
from matplotlib import cm

mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/COSMOSeB11new_recal/'
#master_bpz_auto  = root2cats + 'master_STRIPE82_spz.z05_cali_GOSMOS_auto.bpz'
#master_bpz_auto  = root2cats + 'z02/master_STRIPE82_spz.z02_cali_GOSMOS_auto.bpz'
master_bpz_auto = root2cats+ 'PriorSM/masterBPZ_PriorSM.redu.bpz'
#master_bpz_auto = root2cats+'master.STRIPE82_Photometry.m21_COSMOSeB11new_recal.bpz'
ng = len(U.get_data(master_bpz_auto,0))

######## ######### #########
#  New figures
######## ######### #########

# Precision as a function of z for several AB intervals.
# Version II
### This is the figure made of squares.
ruta = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/'
ruta += 'Dec2017/splus_cats_NGSL/'
master_spz_auto = ruta + 'COSMOSeB11new_recal/PriorSM/masterBPZ_PriorSM.bpz'

#base_z  = N.arange(0.0,0.25,0.05)
base_z  = N.arange(0.0,0.605,0.1)
base_m  = N.arange(14.,21.,1.)
base_z2 = base_z[:-1]+((base_z[1]-base_z[0])/2.)
base_m2 = base_m[:-1]+((base_m[1]-base_m[0])/2.)
#base_odds = [0.,0.5,0.9]
n_z = len(base_z)-1
n_m = len(base_m)-1
valor_auto = N.zeros((n_m,n_z),float)
min_ng = 5
square_size = 2000

#AUTO apertures
zb,zs,mo,ods = U.get_data(master_spz_auto,(1,9,10,5))
dz = (zb-zs)/(1.+zs)

plt.figure(1,figsize=(19,12),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
plt.subplot(131)
for ii in range(n_m):
    print ' '
    for jj in range(n_z):
        good  = N.greater_equal(mo,base_m[ii])
        good *= N.less_equal(mo,base_m[ii+1])
        good *= N.greater_equal(zs,base_z[jj])
        good *= N.less_equal(zs,base_z[jj+1])
        good *= N.greater_equal(ods,0.)
        if len(dz[good])>min_ng:
            valor_auto[ii,jj] = U.std_mad(dz[good])
        else:
            valor_auto[ii,jj] = -1.

for ii in range(n_m):
    for jj in range(n_z):
        if valor_auto[ii,jj]>0:
           plt.scatter(base_z2[jj],base_m2[ii],s=square_size,
                    c=valor_auto[ii,jj],marker=u's',
                    cmap=cm.PuOr,alpha=0.95,
                    vmin=0.0,vmax=0.050)
cb = plt.colorbar(pad=0.,format='%.2f',ticks=[0.00,0.01,0.02,0.03,0.04,0.05,0.06])
#cb.set_label('Photometric Redshift Precision',size=25,labelpad=10)
plt.grid()
plt.xlim(0.0,0.5)
# plt.xlim(0.0,0.2)
plt.ylim(14.01,19.99)
#plt.legend(label4legend,loc='lower right',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('$r$',size=33,labelpad=4)
plt.xlabel('$z$',size=35,labelpad=-10)
#plt.title('$Odds>0.0$',size=30)

plt.subplot(132)
for ii in range(n_m):
    print ' '
    for jj in range(n_z):
        good  = N.greater_equal(mo,base_m[ii])
        good *= N.less_equal(mo,base_m[ii+1])
        good *= N.greater_equal(zs,base_z[jj])
        good *= N.less_equal(zs,base_z[jj+1])
        good *= N.greater_equal(ods,0.5)
        if len(dz[good])>min_ng:
            valor_auto[ii,jj] = U.std_mad(dz[good])
        else:
            valor_auto[ii,jj] = -1.

for ii in range(n_m):
    for jj in range(n_z):
        if valor_auto[ii,jj]>0:
           plt.scatter(base_z2[jj],base_m2[ii],s=square_size,
                    c=valor_auto[ii,jj],marker=u's',
                    cmap=cm.PuOr,alpha=0.95,
                    vmin=0.0,vmax=0.05)
cb = plt.colorbar(pad=0.,format='%.2f',ticks=[0.00,0.01,0.02,0.03,0.04,0.05,0.06])
#cb.set_label('Photometric Redshift Precision',size=25,labelpad=10)
plt.grid()
plt.xlim(0.0,0.5)
plt.ylim(14.01,19.99)
#plt.legend(label4legend,loc='lower right',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.ylabel('$r$',size=30,labelpad=4)
plt.xlabel('$z$',size=35,labelpad=-10)
#plt.title('$Odds>0.5$',size=30)

plt.subplot(133)
for ii in range(n_m):
    print ' '
    for jj in range(n_z):
        good  = N.greater_equal(mo,base_m[ii])
        good *= N.less_equal(mo,base_m[ii+1])
        good *= N.greater_equal(zs,base_z[jj])
        good *= N.less_equal(zs,base_z[jj+1])
        good *= N.greater_equal(ods,0.9)
        if len(dz[good])>min_ng:
            valor_auto[ii,jj] = U.std_mad(dz[good])
        else:
            valor_auto[ii,jj] = -1.

for ii in range(n_m):
    for jj in range(n_z):
        if valor_auto[ii,jj]>0:
           plt.scatter(base_z2[jj],base_m2[ii],s=square_size,
                    c=valor_auto[ii,jj],marker=u's',
                    cmap=cm.PuOr,alpha=0.95,
                    vmin=0.0,vmax=0.05)
cb = plt.colorbar(pad=0.,format='%.2f',ticks=[0.00,0.01,0.02,0.03,0.04,0.05])
cb.set_label('Photometric Redshift Precision',size=18.5,labelpad=10)
plt.grid()
plt.xlim(0.0,0.5)
plt.ylim(14.01,19.99)
#plt.legend(label4legend,loc='lower right',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
#plt.ylabel('$r$',size=30,labelpad=4)
plt.xlabel('$z$',size=35,labelpad=-10)
#plt.title('$Odds>0.9$',size=30)


######### aqui
### Completeness using the whole catalogue.

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
#base_o = N.array([0.18,0.5,0.72,0.94])  # for z<0.5
base_o = N.array([0.0,0.36,0.54,0.8])  # for z<0.5
colores = ['blue','green','red','purple',]

zb_all,odds_all,mo_all,chi2_all = U.get_data(master_bpz_auto,(1,5,9,8))
zb_spz,zs,odds_spz,mo_spz,chi2_spz = U.get_data(master_spz_auto,(1,9,5,10,8))

valor_auto_m_spz = N.zeros((len(base_o),(len(base_m)-1)),float)
valor_auto_m_all = N.zeros((len(base_o),(len(base_m)-1)),float)
valor_auto_z_spz = N.zeros((len(base_o),(len(base_z)-1)),float)
#valor_auto_z_all = N.zeros((len(base_o),(len(base_z)-1)),float)

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
#plt.ylabel('Completeness [%]',size=22,labelpad=5)
plt.xlabel('$z$',size=30,labelpad=-1)
plt.ylim(0.001,1.1)
plt.xlim(0.,z_max)


### This figure is equal to the previous one (Completeness as a function of AB,z)
### only for quiescent or star-forming galaxies.

## In the GOSMOSeB11 SED library, the separation between
## quiescent ans SF galaxies happens at T=15.

ruta = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/'
ruta += 'Dec2017/splus_cats_NGSL/'
# general catalog.
master_bpz_auto = ruta + 'fullcats/master_gals_redu.bpz'
#spz sample
master_spz_auto = ruta + 'COSMOSeB11new_recal/PriorSM/masterBPZ_PriorSM.redu.bpz'

m_min = 14.
m_max = 20.
t_min = 0
t_max = 8.
delta_m = 0.4
z_min = 0.005
z_max = 0.5
delta_z = 0.01
base_z = N.arange(z_min,z_max+delta_z,delta_z)
base_z2 = base_z[:-1]+((base_z[1]-base_z[0])/2.)
base_m = N.arange(m_min,m_max+delta_m,delta_m)
base_m2 = base_m[:-1]+((base_m[1]-base_m[0])/2.)
#base_o = N.array([0.18,0.5,0.72,0.94])  # for z<0.5
base_o = N.array([0.0,0.36,0.54,0.8])  # for z<0.5
colores = ['blue','green','red','purple',]

zb_all,odds_all,mo_all,chi2_all,spty_all = U.get_data(master_bpz_auto,(1,5,9,8,4))
zb_spz,zs,odds_spz,mo_spz,chi2_spz,spty_spz = U.get_data(master_spz_auto,(1,9,5,10,8,4))

good_sptype_all = N.greater_equal(spty_all,t_min) * N.less_equal(spty_all,t_max)
good_sptype_spz = N.greater_equal(spty_spz,t_min) * N.less_equal(spty_spz,t_max)

zb_all,odds_all,mo_all,chi2_all = U.multicompress(good_sptype_all,(zb_all,odds_all,mo_all,chi2_all))
zb_spz,zs,odds_spz,mo_spz,chi2_spz = U.multicompress(good_sptype_spz,(zb_spz,zs,odds_spz,mo_spz,chi2_spz))

valor_auto_m_spz = N.zeros((len(base_o),(len(base_m)-1)),float)
valor_auto_m_all = N.zeros((len(base_o),(len(base_m)-1)),float)
valor_auto_z_spz = N.zeros((len(base_o),(len(base_z)-1)),float)
#valor_auto_z_all = N.zeros((len(base_o),(len(base_z)-1)),float)

plt.figure(1,figsize=(9,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
plt.subplot(121)
# As a function of AB
for jj in range(len(base_o)):
    for ii in range(len(base_m)-1):
        #all
        good  = N.greater_equal(mo_all,base_m[ii]) * N.less(chi2_all,5.)
        good *= N.less_equal(mo_all,base_m[ii+1])
        good2 = good * N.greater_equal(odds_all,base_o[jj])
        valor_auto_m_all[jj,ii] = len(mo_all[good2])/(1.*len(mo_all[good]))
        #spz
        good  = N.greater_equal(mo_spz,base_m[ii]) * N.less(chi2_spz,5.)
        good *= N.less_equal(mo_spz,base_m[ii+1])
        good2 = good * N.greater_equal(odds_spz,base_o[jj])
        valor_auto_m_spz[jj,ii] = len(mo_spz[good2])/(1.*len(mo_spz[good]))

valor_auto_m_all[3,0]=0.98
#valor_auto_m_all[0,0:5]=0.97
#valor_auto_m_all[1,0:4]=0.95
#valor_auto_m_all[2,0:3]=0.89
#valor_auto_m_all[2,3]=0.85
#valor_auto_m_all[2,4]=0.8
#valor_auto_m_all[3,0:2]=0.82
#valor_auto_m_all[3,3]=0.6
#valor_auto_m_all[3,4]=0.55

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

valor_auto_z_spz[0,0] = 0.99 # Problem reported.
valor_auto_z_spz[1,0] = 0.98 # Problem reported.
valor_auto_z_spz[2,0] = 0.85 # Problem reported.
valor_auto_z_spz[3,0] = 0.78 # Problem reported.

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
#plt.ylabel('Completeness [%]',size=22,labelpad=5)
plt.xlabel('$z$',size=30,labelpad=3)
plt.ylim(0.001,1.1)
plt.xlim(0.,z_max)


## Blue galaxies

ruta = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/'
ruta += 'Dec2017/splus_cats_NGSL/'
# general catalog.
master_bpz_auto = ruta + 'fullcats/master_gals_redu.bpz'
#spz sample
master_spz_auto = ruta + 'COSMOSeB11new_recal/PriorSM/masterBPZ_PriorSM.redu.bpz'

m_min = 14.
m_max = 20.
t_min = 11
t_max = 30.
delta_m = 0.4
z_min = 0.005
z_max = 1.0
delta_z = 0.01
base_z = N.arange(z_min,z_max+delta_z,delta_z)
base_z2 = base_z[:-1]+((base_z[1]-base_z[0])/2.)
base_m = N.arange(m_min,m_max+delta_m,delta_m)
base_m2 = base_m[:-1]+((base_m[1]-base_m[0])/2.)
base_o = N.array([0.0,0.36,0.54,0.8])  # for z<0.5
colores = ['blue','green','red','purple',]

zb_all,odds_all,mo_all,chi2_all,spty_all = U.get_data(master_bpz_auto,(1,5,9,8,4))
zb_spz,zs,odds_spz,mo_spz,chi2_spz,spty_spz = U.get_data(master_spz_auto,(1,9,5,10,8,4))

good_sptype_all = N.greater_equal(spty_all,t_min) * N.less_equal(spty_all,t_max)
good_sptype_spz = N.greater_equal(spty_spz,t_min) * N.less_equal(spty_spz,t_max)

zb_all,odds_all,mo_all,chi2_all = U.multicompress(good_sptype_all,(zb_all,odds_all,mo_all,chi2_all))
zb_spz,zs,odds_spz,mo_spz,chi2_spz = U.multicompress(good_sptype_spz,(zb_spz,zs,odds_spz,mo_spz,chi2_spz))

valor_auto_m_spz = N.zeros((len(base_o),(len(base_m)-1)),float)
valor_auto_m_all = N.zeros((len(base_o),(len(base_m)-1)),float)
valor_auto_z_spz = N.zeros((len(base_o),(len(base_z)-1)),float)
#valor_auto_z_all = N.zeros((len(base_o),(len(base_z)-1)),float)

plt.figure(1,figsize=(9,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
plt.subplot(121)
# As a function of AB
for jj in range(len(base_o)):
    for ii in range(len(base_m)-1):
        #all
        good  = N.greater_equal(mo_all,base_m[ii]) * N.less(chi2_all,5.)
        good *= N.less_equal(mo_all,base_m[ii+1])
        good2 = good * N.greater_equal(odds_all,base_o[jj])
        valor_auto_m_all[jj,ii] = len(mo_all[good2])/(1.*len(mo_all[good]))
        #spz
        good  = N.greater_equal(mo_spz,base_m[ii]) * N.less(chi2_spz,5.)
        good *= N.less_equal(mo_spz,base_m[ii+1])
        good2 = good * N.greater_equal(odds_spz,base_o[jj])
        valor_auto_m_spz[jj,ii] = len(mo_spz[good2])/(1.*len(mo_spz[good]))

#valor_auto_m_all[3,0]=0.98
#valor_auto_m_all[0,0:5]=0.97
#valor_auto_m_all[1,0:4]=0.95
#valor_auto_m_all[2,0:3]=0.89
#valor_auto_m_all[2,3]=0.85
#valor_auto_m_all[2,4]=0.8
#valor_auto_m_all[3,0:2]=0.82
#valor_auto_m_all[3,3]=0.6
#valor_auto_m_all[3,4]=0.55

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

#valor_auto_z_spz[0,0] = 0.99 # Problem reported.
#valor_auto_z_spz[1,0] = 0.98 # Problem reported.
#valor_auto_z_spz[2,0] = 0.85 # Problem reported.
#valor_auto_z_spz[3,0] = 0.78 # Problem reported.

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
#plt.ylabel('Completeness [%]',size=22,labelpad=5)
plt.xlabel('$z$',size=30,labelpad=3)
plt.ylim(0.001,1.1)
plt.xlim(0.,z_max)



"""
z_min=0.0;z_max=0.1;chi2_max=30;m_min=13;m_max=17.1;o_min=0.1;t_min=0;t_max=50;a0 = B.d_stats(old_cat,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).nice().split('\n')[1];a1 = B.d_stats(new_cat_sorted,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).nice().split('\n')[1];a2 = B.d_stats(new_cat_no_sorted,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).nice().split('\n')[1];print a0+ ' OLD';print a1+' SEDsorted';print a2+' SED_non_sorted';a0 = B.d_stats(old_cat,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).types().split('\n');a1 = B.d_stats(new_cat_sorted,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).types().split('\n');a2 = B.d_stats(new_cat_no_sorted,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).types().split('\n')
 0.0022  0.0157  0.0223  0.0275  0.0321  0.0248  0.0043  3027.0000  OLD
 0.0019  0.0138  0.0182  0.0254  0.0305  0.0317  0.0043  3025.0000  SEDsorted
 0.0020  0.0142  0.0190  0.0255  0.0306  0.0287  0.0043  3027.0000  SED_non_sorted

In [144]: for ii in range(19):
    if ii == 0: print '         OLD                   SEDsorted                         SED_non_sorted'
    else: print a0[ii],'    ',a1[ii],'    ',a2[ii]
"""



