__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import matplotlib.pyplot as plt
import numpy as N

# Root to data
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/COSMOSeB11new_recal/'
bpz_cat = root2cats + 'master.STRIPE82_Photometry.m21_COSMOSeB11new_recal_redu.bpz'

# Reading data
zb,mo,od,zs,chi2,tb = U.get_data(bpz_cat,(1,10,5,9,8,4))

# Doing some cleaning
good = N.less(chi2,30) * N.greater_equal(mo,14)
good *= N.greater_equal(od,0.05) * N.greater_equal(zs,0.001)
zb,mo,od,zs,chi2,tb = U.multicompress(good,(zb,mo,od,zs,chi2,tb))
dz = (zs-zb)/(1.+zs)

#base magnitudes
basem  = N.arange(14,22,1.)
basem2 = basem[:-1]+((basem[1]-basem[0])/2.)

# base templates
baset  = N.arange(1,max(tb)+1,3)
baset2 = baset[:-1]+((baset[1]-baset[0])/2.)

# base redshift
basez  = N.arange(0.,0.6,0.1)
basez2 = basez[:-1]+((basez[1]-basez[0])/2.)

#Odds base
odds = N.array([0.18,0.5,0.72,0.94])
odds2 = N.array([0.18,0.5,0.72,0.94,1.00])
# colors
colores = ['blue','green','red','purple',]

#plots
res = 20
mor = mo[::res]
dzr = dz[::res]
odr = od[::res]
tbr = tb[::res]
zsr = zs[::res]
dzr = dz[::res]

plt.figure(10,figsize=(8,12),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
for ss in range(18):
    plt.subplot(6,3,ss+1)
    if ss==0:
       plt.ylabel('$\Delta_{z}/(1+z_{s})$',size=15)
       plt.xticks([])
       #plt.yticks([])
       for ii in range(4):
           linea_base = U.bin_stats(mor[odr>odds[ii]],dzr[odr>odds[ii]],basem)
           plt.plot(basem,linea_base,'-s',color=colores[ii],ms=6,alpha=0.2)
       plt.ylim(-0.05,0.05)
    if ss==1:
       #plt.ylabel('$\Delta_{z}/(1+z_{s})$',size=15)
       plt.xticks([])
       plt.yticks([])
       for ii in range(4):
           linea_base = U.bin_stats(tbr[odr>odds[ii]],dzr[odr>odds[ii]],baset)
           plt.plot(baset,linea_base,'-s',color=colores[ii],ms=6,alpha=0.2)
       plt.ylim(-0.05,0.05)
    if ss==2:
       plt.xticks([])
       plt.yticks([])
       for ii in range(4):
           linea_base = U.bin_stats(zsr[odr>odds[ii]],dzr[odr>odds[ii]],basez)
           plt.plot(basez,linea_base,'-s',color=colores[ii],ms=6,alpha=0.2)
       plt.ylim(-0.05,0.05)

    if ss==3:
       plt.ylabel('$\sigma_{z}$ $x10^{-2}$',size=22)
       plt.xticks([])
       #plt.yticks([])
       for ii in range(4):
           sigma_zm,bb = get_sigma_z(dzr[odr>odds[ii]],mor[odr>odds[ii]],basem)
           plt.plot(bb,sigma_zm*100.,'-s',color=colores[ii],ms=6,alpha=0.2)
       plt.ylim(0.0,4.)

    if ss==4:
       plt.xticks([])
       plt.yticks([])
       for ii in range(4):
           sigma_zt,bb = get_sigma_z(dzr[odr>odds[ii]],tbr[odr>odds[ii]],baset)
           plt.plot(bb,sigma_zt*100.,'-s',color=colores[ii],ms=6,alpha=0.2)
       plt.ylim(0.,4.)

    if ss==5:
       plt.xticks([])
       plt.yticks([])
       for ii in range(4):
           sigma_zs,bb = get_sigma_z(dzr[odr>odds[ii]],zsr[odr>odds[ii]],basez)
           plt.plot(bb,sigma_zs*100.,'-s',color=colores[ii],ms=6,alpha=0.2)
       plt.ylim(0.,4.)

    if ss==6:
       plt.ylabel('Bias $x10^{-2}$',size=15)
       plt.xticks([])
       for ii in range(4):
           sigma_zm,bb = get_mean_z(dzr[odr>odds[ii]],mor[odr>odds[ii]],basem)
           plt.plot(bb,sigma_zm*100.,'-s',color=colores[ii],ms=6,alpha=0.2)
       plt.ylim(-0.5,0.5)

    if ss==7:
       plt.xticks([])
       plt.yticks([])
       for ii in range(4):
           sigma_zt,bb = get_mean_z(dzr[odr>odds[ii]],tbr[odr>odds[ii]],baset)
           plt.plot(bb,sigma_zt*100.,'-s',color=colores[ii],ms=6,alpha=0.2)
       plt.ylim(-0.5,0.5)

    if ss==8:
       plt.xticks([])
       plt.yticks([])
       for ii in range(4):
           sigma_zs,bb = get_mean_z(dzr[odr>odds[ii]],zsr[odr>odds[ii]],basez)
           plt.plot(bb,sigma_zs*100,'-s',color=colores[ii],ms=6,alpha=0.2)
       plt.ylim(-0.5,0.5)

    if ss==9:
       plt.ylabel('Compl',size=15)
       plt.xticks([])
       mean_zs,bb = get_completeness(dzr,zsr,basez,odr,odds2)
       #mean_zs,bb = get_completeness(dzr[odr>odds[ii]],zsr[odr>odds[ii]],basez)
       for ii in range(4):
           plt.plot(bb,mean_zs[:,ii]*100,'-s',color=colores[ii],ms=6,alpha=0.2)

    if ss==15:
       plt.xlabel('$r$',size=22)
       plt.xticks([])

    if ss==16:
       plt.xlabel('$t$',size=22)
       plt.xticks([])
       plt.yticks([])


    if ss==17:
       plt.xlabel('$z$',size=22)
       plt.xticks([])
       plt.yticks([])






def get_sigma_z(ddz,vector,base):

    n_ele = len(base)-1
    values = N.zeros(n_ele)
    base2 = base[:-1]+((base[1]-base[0])/2.)
    for ii in range(n_ele):
        good = N.greater_equal(vector,base[ii])
        good *= N.less_equal(vector,base[ii+1])
        values[ii] = U.std_mad(ddz[good])
    return values,base2


def get_mean_z(ddz,vector,base):

    n_ele = len(base)-1
    values = N.zeros(n_ele)
    base2 = base[:-1]+((base[1]-base[0])/2.)
    for ii in range(n_ele):
        good = N.greater_equal(vector,base[ii])
        good *= N.less_equal(vector,base[ii+1])
        values[ii] = U.mean_robust(ddz[good])
    return values,base2

def get_completeness(ddz,vector,base,vectorod,baso):
    n_ele = len(base)-1
    n_ods = len(baso)-1
    values = N.zeros((n_ele,n_ods),float)
    base2 = base[:-1]+((base[1]-base[0])/2.)
    for jj in range(n_ods):
        good_ods = N.greater_equal(vectorod,baso[jj])
        good_ods *= N.less_equal(vectorod,baso[jj+1])
        vector_redu = vector[good_ods]
        ddzr = ddz[good_ods]
        for ii in range(n_ele):
            good = N.greater_equal(vector_redu,base[ii])
            good *= N.less_equal(vector_redu,base[ii+1])
            values[ii,jj] = len(ddzr[good])/(len(ddzr)*1.)

    return values,base2

