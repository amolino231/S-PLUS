__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt

root_to_cat = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/'
root_to_cat += 'Dec2017/splus_cats_NGSL/'
photo_cat = root_to_cat+'cat/master.STRIPE82_Photometry.m21.cat'
u,g,r,i,z,zs = U.get_data(photo_cat,(15,60,78,96,114,125))
bpz_6 = root_to_cat + 'COSMOSeB11new_recal/master.STRIPE82_Photometry.m21_COSMOSeB11new_recal.bpz'
zb,mo,od,zs,chi2,tb = U.get_data(bpz_6,(1,10,5,9,8,4))
dz = (zb-zs)/(1.+zs)
red = N.less_equal(tb,2.5) * N.greater(zs,0.0) * N.less_equal(zs,0.3) * N.less(abs(r),19.)
z_threshold = 0.03
good_Red = red * N.less_equal(abs(dz), z_threshold)
bad_Red = red * N.greater(abs(dz), z_threshold)
color_x = u-g
color_y = g-r
color_x2 = r-i
color_y2 = i-z
dd=0.05
plt.clf()
plt.subplot(221)
a1,a2,a3 = plt.hist(color_x[bad_Red],N.arange(-0.5,3.0,2.*dd),facecolor='red',alpha=0.2)
a1,a2,a3 = plt.hist(color_x[good_Red],N.arange(-0.5,3.0,2.*dd),facecolor='blue',alpha=0.2)
plt.xlim(0.,3.)
plt.grid()
plt.xlabel("u-g",size=20)
plt.legend(['bad','good'],loc='best',numpoints=1,fontsize=15)
plt.subplot(222)
a1,a2,a3 = plt.hist(color_x2[bad_Red],N.arange(-0.5,3.0,dd),facecolor='red',alpha=0.2)
a1,a2,a3 = plt.hist(color_x2[good_Red],N.arange(-0.5,3.0,dd),facecolor='blue',alpha=0.2)
plt.xlabel("g-r",size=20)
plt.grid()
plt.xlim(-0.2,1.)
plt.legend(['bad','good'],loc='best',numpoints=1,fontsize=15)
plt.subplot(223)
a1,a2,a3 = plt.hist(color_x2[bad_Red],N.arange(-0.5,3.0,dd),facecolor='red',alpha=0.2)
a1,a2,a3 = plt.hist(color_x2[good_Red],N.arange(-0.5,3.0,dd),facecolor='blue',alpha=0.2)
plt.xlabel("r-i",size=20)
plt.grid()
plt.xlim(-0.2,1.)
plt.legend(['bad','good'],loc='best',numpoints=1,fontsize=15)
plt.subplot(224)
a1,a2,a3 = plt.hist(color_y2[bad_Red],N.arange(-0.5,3.0,dd),facecolor='red',alpha=0.2)
a1,a2,a3 = plt.hist(color_y2[good_Red],N.arange(-0.5,3.0,dd),facecolor='blue',alpha=0.2)
plt.xlabel("i-z",size=20)
plt.xlim(-0.3,0.5)
plt.grid()
plt.legend(['bad','good'],loc='best',numpoints=1,fontsize=15)


