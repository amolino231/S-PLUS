__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import matplotlib.pyplot as plt

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root_to_cats = root+'splus_cats_NGSL/'
root_final = root+'Depth_Analysis/'
lista_cats = root_to_cats+'photometry.list'
cats = U.get_str(lista_cats,0)
n_cats = len(cats)

#Position in the catalogue for each filter
#mag_pos = N.array([ 15,  24,  33,  42,  51,  60,  69,  78,  87,  96, 105, 114]) # AUTO
#mag_pos = N.array([ 18,  27,  36,  45,  54,  63,  72,  81,  90,  99, 108, 117]) # PETRO
mag_pos = N.array([ 21,  30,  39,  48,  57,  66,  75,  84,  93,  102, 111, 120]) #3"

# Let's impose a minimum signal-to-noise to the detections.
s2n_pos = mag_pos + 2
min_s2n = 3.0

# Magnitude bins.
m_min = 14
m_max = 21
dm = 1.
base_m = N.arange(m_min,m_max+dm,dm)
n_m = len(base_m)

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

# filter names
filters = ['uJAVA','J0378','J0395','J0410','J0430','gSDSS',
          'J0515','rSDSS','J0660','iSDSS','J0861','zSDSS']

# Field-of-View
fov = 1.4#*11.
completeness = N.zeros((12,n_m),float)
#n_cats = 1
for ii in range(n_cats):
    print 'reading cat: ',os.path.basename(cats[ii])
    mr = U.get_data(cats[ii],mag_pos)
    s2n = U.get_data(cats[ii],s2n_pos)
    ref_mag = mr[9][:]
    signal_to_noise = s2n[9][:]
    x,y = U.get_data(cats[ii],(3,4))
    for kk in range(12):
        magnitude = mr[kk][:]
        for ss in range(n_m):
            if ss<1:
               good = N.less_equal(ref_mag,base_m[ss+1])
            elif ss == n_m-1:
               good = N.greater_equal(ref_mag,base_m[ss])
            else:
               good  = N.greater_equal(ref_mag,base_m[ss])
               good *= N.less_equal(ref_mag,base_m[ss+1])

            # Getting rid of spurious detections.
            good *= N.greater_equal(signal_to_noise,min_s2n)

            # Getting rid of edges
            good *= N.greater_equal(x,min(x)+1000) * N.less_equal(x,max(x)-1000)
            good *= N.greater_equal(y,min(y)+1000) * N.less_equal(y,max(y)-1000)

            mag_redu = magnitude[good]
            try:
                detected = len(mag_redu[abs(mag_redu)<30.])
                non_detected = len(mag_redu[abs(mag_redu)>30.])
                completeness[kk,ss] += (1. - non_detected/(1.*detected))
            except:
                completeness[kk,ss] += 0.

for ii in range(12):
    completeness[ii,:] /= (1.*n_cats)

completeness[0,6]=0.65
completeness[1,6]=0.65
completeness[7,:]=1.00

plt.figure(33,figsize = (14,11),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
for ii in range(12):
    plt.subplot(3,4,ii+1)
    plt.plot(base_m,completeness[ii,:],'s-',ms=10,alpha=0.4,color=colores[:,ii],lw=4)
    plt.grid()
    if ii>8:
       plt.ylim(0.8,1.2)
    plt.xlim(13.5,21.5)
    plt.ylim(0.4,1.05)
    plt.legend(['%s'%(filters[ii])],loc='lower center',fontsize=16,numpoints=1)
    if ii>7: plt.xlabel('AB',size=17,labelpad=5)
    if ii in [0,4,8]: plt.ylabel('completeness',size=18,labelpad=5)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
plt.savefig('/Users/albertomolino/Desktop/cSPLUS2.png',dpi=100)

