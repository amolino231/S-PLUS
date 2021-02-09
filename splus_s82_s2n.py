__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import bpz_tools as B
import matplotlib.pyplot as plt

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root_to_cats = root + 'splus_cats_NGSL/'
lista_cats = root_to_cats+'photometry.list'
cats = U.get_str(lista_cats,0)
n_cats = len(cats)

#Position in the catalogue for each filter
#mag_pos = N.array([ 21,  30,  39,  48,  57,  66,  75,  84,  93,  102, 111, 120])
mag_pos = N.array([ 18,  27,  36,  45,  54,  63,  72,  81,  90,  99, 108, 117])
s2n_pos = mag_pos+2
#s2n_pos = N.array([ 20,  29,  38,  47,  56,  65,  74,  83,  92,  101, 110, 119])

#colors for different filters
colores = N.zeros((3,12),float)
colores[:,0]=(0.00,0.00,1.00)
colores[:,1]=(0.00,0.25,1.00)
colores[:,2]=(0.00,0.65,1.00)
colores[:,3]=(0.00,0.50,0.00)
colores[:,4]=(0.85,0.65,0.00)
colores[:,5]=(0.75,0.50,0.00)
colores[:,6]=(0.80,0.25,0.00)
colores[:,7]=(1.00,0.00,0.00)
colores[:,8]=(0.85,0.00,0.00)
colores[:,9]=(0.65,0.00,0.00)
colores[:,10]=(0.35,0.00,0.00)
colores[:,11]=(0.25,0.00,0.00)

#basem
basem = N.arange(14,21.5,0.5)
nm = len(basem)

#resutls
results = N.zeros((nm,12),float)

# filter names
filters = ['uJAVA','J0378','J0395','J0410','J0430','gSDSS',
          'J0515','rSDSS','J0660','iSDSS','J0861','zSDSS']
for ii in range(n_cats):
    print 'reading catalogue: ',cats[ii]
    mags = U.get_data(cats[ii],mag_pos)
    s2ns = U.get_data(cats[ii],s2n_pos)
    for kk in range(12):
        temporal_mag = mags[kk][:]
        temporal_s2n = s2ns[kk][:]
        results[:,kk] += B.bin_stats(temporal_mag,temporal_s2n,basem,'mean_robust')


results[-1,0]=174
results[-1,1]=136
results[-1,2]=80
results[-1,3]=80
results[-1,4]=90
results[-1,6]=110
results[-1,8]=300
results[-1,9]=250
results[-1,10]=95
results[-1,11]=130


plt.figure(10,figsize = (14,11),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
for iii in range(12):
    plt.subplot(3,4,iii+1)
    plt.semilogy(basem,results[:,iii]/(1.*n_cats),'-',lw=12,color=colores[:,iii],alpha=0.5)
    plt.legend([filters[iii]],fontsize=17,loc='lower left')
    plt.semilogy(basem,results[:,5]/(1.*n_cats),'--',color='black',lw=6,alpha=0.5)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if iii>7: plt.xlabel('AB',size=17,labelpad=5)
    if iii in [0,4,8]: plt.ylabel('signal-to-noise',size=18,labelpad=3)
    plt.xlim(13,23)
    plt.ylim(0.1,2000)
    #plt.ylim(10.,max(a10/fov)*1.2)
    #plt.legend(['$S/N$ $>1:$    $%.2f$'%(max0),'$S/N$ $>3:$    $%.2f$'%(max1),'$S/N$ $>5:$    $%.2f$'%(max2),'$S/N$ $>10:$   $%.2f$'%(max3)],loc='upper left',fontsize=14)
    # plt.legend(['$S/N$ $>0:$    $%.1f$'%(max0),'$S/N$ $>3:$    $%.1f$'%(max1),'$S/N$ $>5:$    $%.1f$'%(max2),'$S/N$ $>10:$   $%.1f$'%(max3),'$S/N$ $>100:$ $%.1f$'%(max4),],loc='upper left',fontsize=25)
    plt.grid()
