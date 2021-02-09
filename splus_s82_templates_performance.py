__author__ = 'albertomolino'

"""
This routine creates a mosaic-like plot
showing per each SED model the observed
precision and bias using the S82 data.
---
It needs all the individual BPZ files.

"""

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import bpz_tools as B
import useful as U
import numpy as N
import matplotlib.pyplot as plt

bpz_path = '/Users/albertomolino/codigos/bpz-1.99.2/'
sed_path = bpz_path+'SED/'
sed_lib = 'COSMOSeB11new_recal.list'
sed_models = U.get_str(sed_path+sed_lib,0)
n_models = len(sed_models)
if n_models <> 12: n_models = 12

bpzlist = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/splus_cats_NGSL/'
#bpzlist += 'COSMOSeB11new/master.STRIPE82_Photometry.m21_COSMOSeB11new.bpz'
#bpzlist += 'COSMOSeB11new_recal/bpz_short.list'
bpzlist += 'COSMOSeB11new_recal/PriorSM/masterBPZ_PriorSM.list'
bpz_files  = U.get_str(bpzlist,0)
n_bpzs = len(bpz_files)

## 1 accuracy and bias of each template.

med_value = N.zeros((n_bpzs,n_models),float)
std_value = N.zeros((n_bpzs,n_models),float)
out_value = N.zeros((n_bpzs,n_models),float)
num_sourc = N.zeros((n_bpzs,n_models),float)
for ii in range(n_bpzs):
    #print bpz_files[ii]
    ao = B.d_stats(bpz_files[ii],mmin=14,mmax=19,zmax=0.5)
    pepe = ao.types().split('\n')[1:-1]
    for jj in range(n_models):
        try:
            pepa = pepe[jj].split()
            #print pepa
            med_value[ii,jj] = float(pepa[1])
            std_value[ii,jj] = float(pepa[2])
            out_value[ii,jj] = float(pepa[3])
            num_sourc[ii,jj] = float(pepa[4])
        except:
            med_value[ii,jj] = 0.0
            std_value[ii,jj] = 0.0
            out_value[ii,jj] = 0.0
            num_sourc[ii,jj] = 0.0

print 'Model  med  std  out  num '
for jj in range(n_models):
    good = N.greater(std_value[:,jj],0.0001)
    linea  = '%i, %.3f, %.3f, '%(jj+1, U.mean_robust(med_value[good,jj]),U.mean_robust(std_value[good,jj]))
    linea += '%.3f  %i '%(U.mean_robust(out_value[good,jj]),U.sum(num_sourc[good,jj]))
    print linea


base_sz = N.arange(0.0005,0.095,0.015)
base_sz_2 = N.arange(-0.095,0.095,0.015)
plt.clf()
for ss in range(n_models):
    plt.subplot(3,4,ss+1)
    good = N.greater(std_value[:,ss],0.0001)
    a1,a2,a3 = plt.hist(std_value[good,ss],base_sz,alpha=0.5,normed=1)
    b1,b2,b3 = plt.hist(med_value[good,ss],base_sz_2,facecolor='red',alpha=0.5,normed=1)
    plt.legend(['$<\sigma>$:%.1f'%(U.mean_robust(std_value[good,ss])*100.),'$<\mu>$:%.1f'%(U.mean_robust(med_value[good,ss])*100.)]
    ,loc='upper left',fontsize=12)
    plt.grid()
    plt.xlim(-0.139,0.099)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    if ss in [8,9,10,11]:
       plt.xlabel('$\Delta_{z}/(1+z)$',size=18,labelpad=12)
    if ss in [0,4,8]:
       plt.ylabel('$\#$',size=18,labelpad=12)


## 2 Study of the symmetry as a function of types, z and mag.
master_bpz_root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82'
master_bpz_root += '/Dec2017/splus_cats_NGSL/COSMOSeB11new_recal/'
master_bpz = master_bpz_root + 'master.STRIPE82_Photometry.m21_COSMOSeB11new_recal_redu.bpz'
zb,zs,mm,tb,od = U.get_data(master_bpz,(1,9,10,4,5))
dz = (zb-zs)/(1.+zs)
pepe = abs(base_sz)
pos = N.where(pepe==min(pepe))[0][0]
base_sz3 = N.arange(-0.095,0.095,0.0001)
## 2.1 as a function of magnitude.
base_m = N.arange(17,22.,1.)
n_bins_m = len(base_m)-1
value_50_m = N.zeros(n_bins_m)
for hh in range(n_bins_m):
    sample  = N.greater_equal(mm,base_m[hh])
    sample *= N.less_equal(mm,base_m[hh+1])
    a1,a2,a3 = plt.hist(dz[sample],base_sz3,cumulative=1,normed=1)
    pepe = abs(a1-0.5)
    pos = N.where(pepe==min(pepe))[0][0]
    value_50_m[hh] = a1[pos]
    print '%.1f <r< %.1f: v50: %.3f (%i)'%(base_m[hh],base_m[hh+1],value_50_m[hh],len(dz[sample]))

"""
17.0 <r< 18.0: v50: 0.501 (6936)
18.0 <r< 19.0: v50: 0.499 (9984)
19.0 <r< 20.0: v50: 0.499 (11286)
20.0 <r< 21.0: v50: 0.500 (4611)
"""


## 2.2 as a function of redshift.
base_z = N.arange(0.,0.6,0.1)
n_bins_z = len(base_z)-1
value_50_z = N.zeros(n_bins_z)
for hh in range(n_bins_z):
    sample  = N.greater_equal(zs,base_z[hh])
    sample *= N.less_equal(zs,base_z[hh+1])
    a1,a2,a3 = plt.hist(dz[sample],base_sz3,cumulative=1,normed=1)
    pepe = abs(a1-0.5)
    pos = N.where(pepe==min(pepe))[0][0]
    value_50_z[hh] = a1[pos]
    print '%.2f <z< %.2f: v50: %.3f (%i)'%(base_z[hh],base_z[hh+1],value_50_z[hh],len(dz[sample]))

"""
0.00 <z< 0.10: v50: 0.500 (8765)
0.10 <z< 0.20: v50: 0.500 (11927)
0.20 <z< 0.30: v50: 0.499 (5729)
0.30 <z< 0.40: v50: 0.500 (4781)
0.40 <z< 0.50: v50: 0.499 (4654)
"""

## 2.3 as a function of Template.
base_t = N.arange(1,17,1)
n_bins_t = len(base_t)-1
value_50_t = N.zeros(n_bins_t)
for hh in range(n_bins_t):
    sample  = N.greater_equal(tb,base_t[hh])
    sample *= N.less_equal(tb,base_t[hh+1])
    a1,a2,a3 = plt.hist(dz[sample],base_sz3,cumulative=1,normed=1)
    pepe = abs(a1-0.5)
    pos = N.where(pepe==min(pepe))[0][0]
    value_50_t[hh] = a1[pos]
    #print len(dz[sample])
    print 'Tb %i: v50: %.3f (%i)'%(base_t[hh],value_50_t[hh],len(dz[sample]))

"""
Tb 1: v50: 0.500 (632)
Tb 2: v50: 0.502 (823)
Tb 3: v50: 0.500 (955)
Tb 4: v50: 0.501 (1409)
Tb 5: v50: 0.501 (5211)
Tb 6: v50: 0.501 (8449)
Tb 7: v50: 0.500 (2656)
Tb 8: v50: 0.501 (1949)
Tb 9: v50: 0.500 (1197)
Tb 10: v50: 0.500 (2446)
Tb 11: v50: 0.500 (3666)
Tb 12: v50: 0.501 (4672)
Tb 13: v50: 0.500 (1638)
Tb 14: v50: 0.498 (427)
Tb 15: v50: 0.516 (47)
"""