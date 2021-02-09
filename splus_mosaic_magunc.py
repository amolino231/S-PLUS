__author__ = 'albertomolino'

"""
This routine creates a mosaic-like figure showing
the magnitude difference between models and data for
a sample of galaxies in every filter.

"""

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import bpz_tools as B
import matplotlib.pyplot as plt
import phz_plots as P
import alhambra_photools as A

# filter names
filters = ['uJAVA','J0378','J0395','J0410','J0430','gSDSS',
          'J0515','rSDSS','J0660','iSDSS','J0861','zSDSS']

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
columns = root + 'splus_auto.columns'
fluxcomparison = root + 'splus_cats_NGSL/COSMOSeB11new_recal_OT/masterOT_COSMOSeB11new_recal_SM.flux_comparison'

mmin=14.
mmax=18.
basem_high = N.arange(-1.0,1.0,0.01)
basem_med  = N.arange(-3.0,3.0,0.02)
basem_low  = N.arange(-3.0,3.0,0.08)
#basem2_high = basem_high[:-1]+((basem_high[1]-basem_high[0])/2.)
#basem2_low  = basem_low[:-1]+((basem_low[1]-basem_low[0])/2.)

# Reading data from catalogs.
ft1,fob1,efob1 = P.getinfo4_mosaic_uncert_obsvsmod(columns,fluxcomparison)
mag1 = U.get_data(fluxcomparison,1)

#emob = B.e_frac2mag(efob1[ii][:])/float(fob1[ii][:])

plt.figure(1, figsize=(24,12),dpi=80, facecolor='w', edgecolor='k')
plt.clf()

bin_labels = [0,6]
base_x = N.linspace(-3.0,3.0,300)

for ss in range(12):
    plt.subplot(2,6,ss+1)
    if ss==0:
       g_1 = U.greater_equal(mag1,mmin) * U.less_equal(mag1,mmax+3)
    else:
       g_1 = U.greater_equal(mag1,mmin) * U.less_equal(mag1,mmax)
    dm2_1 = B.flux2mag(fob1[ss][g_1])-B.flux2mag(ft1[ss][g_1])
    emob = B.e_frac2mag(efob1[ss][g_1])/fob1[ss][g_1]
    sense_values_1 = N.less(abs(dm2_1),5.)
    sense_values_1 *= N.less(abs(emob),1.)
    dm2_1,emob2 = U.multicompress(sense_values_1,(dm2_1,emob))
    #dm2_1 = N.compress(sense_values_1,dm2_1)

    valor = base_x * 0.
    for ii in range(len(emob2)):
        temporal = A.gaussian(base_x,emob2[ss],0.00,1)
        valor += temporal*temporal
    valor = N.sqrt(valor)

    if ss in [0,1,2,3,4]:
        basem = basem_low
        #basem2 = basem2_low
    elif ss in [5,6]:
        basem = basem_med
    else:
        basem = basem_high
        #basem2 = basem2_high

    v1,v2,v3 = plt.hist(dm2_1,basem,facecolor='grey',alpha=0.3,normed=1)

    pos=N.where(v1==max(v1))[0][0]
    if ss != 5: off=v2[pos+1]
    else: off=v2[pos]

    plt.plot(base_x+off,max(v1)*(valor/valor.max()),'r--',lw=2)
    v1,v2,v3 = plt.hist(dm2_1,basem,histtype='step',color='black',alpha=0.7,linewidth=1.5,normed=1)

    label1 = '%s'%(filters[ss])
    plt.title(label1,fontsize=16.5)
    plt.legend(['Noise'],loc='upper center',fontsize=16)

    if ss>5: plt.xlabel('$\delta_{m}$',size=30,labelpad=5)
    plt.yticks(fontsize=20)
    if ss in bin_labels: plt.ylabel('frequence',size=22,labelpad=10)
    if ss in [0,1,2,3,4]:
        plt.xticks([-1.5,0.0,1.5],['-1.5','0.0','1.5'],fontsize=18)
        plt.xlim(-1.99,1.99)
    elif ss in [5,6]:
        plt.xticks([-0.3,0.0,0.3],['-0.3','0.0','0.3'],fontsize=18)
        plt.xlim(-0.39,0.39)
    else:
        plt.xticks(fontsize=20)
        plt.xlim(-0.29,0.29)
        plt.xticks([-0.2,0.0,0.2],['-0.2','0.0','0.2'],fontsize=18)
    #plt.ylim(0.001,N.max(v1)*1.1)
    plt.ylim(0.001,N.max(v1)*1.2)

    plt.yticks([])
    plt.grid()
    plt.ion()
    plt.show()

# Saving figure
#plt.savefig(fluxcomparison[:-15]+'_mosaic_noise.png',dpi=100)