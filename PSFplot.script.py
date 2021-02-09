#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

import os,sys
import numpy as N
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import splus_calib_tools as sct
import matplotlib.pyplot as plt

#from matplotlib import cm
#import redseq as R

#### script to read and normalize the catalogue to the mean PSF.
root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
final_path = root + 'data_quality/'
lista_cats = root+'splus_cats_NGSL/photometry.list'
cats = U.get_str(lista_cats,0)
n_cats = len(cats)

seeing = N.zeros(n_cats)
x_global = []
y_global = []
s2n_global = []
Flag_global = []
fw_global = []
mr_global = []
for ii in range(n_cats):
    catalog = cats[ii]
    print '%i out of %i'%(ii+1,n_cats)
    x,y,s2n,Flag,fw,mr = U.get_data(catalog,(3,4,6,7,8,84))
    seeing[ii],stars = sct.get_seeing_from_data_pro(fw,mr)
    # It selects only stars
    #x,y,s2n,Flag,fw,mr = U.multicompress(stars,(x,y,s2n,Flag,fw,mr))
    n_sources = len(x)
    for kk in range(n_sources):
        x_global.append(x[kk]-x.min())
        y_global.append(y[kk]-y.min())
        s2n_global.append(s2n[kk])
        Flag_global.append(Flag[kk])
        fw_global.append(fw[kk]/(1.*seeing[ii]))
        #fw_global.append(fw[kk]-seeing[ii])
        mr_global.append(mr[kk])

n_total = len(x_global)
file_out_name = final_path+'psfnew_dividedbyseeing.cat'
file_out = open(file_out_name,'w')
file_out.write('# x y s2n Flag fwhm mr_aper \n')
for ss in range(n_total):
    linea = '%i  %i  %.2f  %i  %.2f  %.2f  \n'%(x_global[ss],y_global[ss],s2n_global[ss],
                                                   Flag_global[ss],fw_global[ss],mr_global[ss])
    file_out.write(linea)
file_out.close()

#########


#psfcat = final_path+'psfnew.cat'
psfcat = final_path+'psfnew_dividedbyseeing.cat'
x,y,s2n,Flag,fw,mr = U.get_data(psfcat,(0,1,2,3,4,5))
good = N.greater(s2n,100)*N.greater(mr,13.5)*N.less(mr,18.5)*N.less(Flag,1)
#good *= N.greater(fw,0.92)*N.less(fw,1.05)
x,y,s2n,Flag,fw,mr = U.multicompress(good,(x,y,s2n,Flag,fw,mr))
rad = N.sqrt((x-N.mean(x))*(x-N.mean(x))+(y-N.mean(y))*(y-N.mean(y)))
rad2=rad*0.55/60.
plt.figure(1)
mat,v1,v2 = R.CC_numberdensity_contour_rangefixed(fw,rad2,0.1,0,1,0)
base=N.arange(0.8,1.2,0.1)

plt.figure(2,figsize = (17,10),dpi=70, facecolor='w', edgecolor='k')
plt.clf()

dos = plt.axes([0.35,0.215,0.2,0.7])
a1,a2,a3 = plt.hist(fw2,base,facecolor='blue',alpha=0.15,orientation='horizontal',normed=1)
a1,a2,a3 = plt.hist(fw2,base,histtype='step',color='blue',alpha=0.5,orientation='horizontal',normed=1)
plt.ylim(0.81,1.19)
plt.grid()
plt.xticks(fontsize=20)
plt.xlabel('frequence',size=25,labelpad=15)
plt.xlim(0,12)
# ;plt.legend(['$\#$4000'],loc='upper right',fontsize=12)

uno = plt.axes([0.1,.215,0.25,0.7])
plt.contour(v1,v2,mat,500,linewidths=2)
plt.xlim(0,60)
plt.ylabel('normed FWHM',size=25,labelpad=15)
plt.xlabel('radial distance [arcmin]',size=25,labelpad=15)
plt.xticks(fontsize=20);plt.yticks(fontsize=20)
plt.ylim(0.81,1.19)
plt.grid()
plt.xlim(1,59)

tres = plt.axes([0.625,0.215,0.3,0.7])
plt.scatter(x/1000., y/1000.,s=120,c=fw, marker=u'o', cmap=cm.gist_rainbow_r, alpha=0.2,vmin=0.9,vmax=1.06)
cb = plt.colorbar(pad=0.,format='%.2f')
cb.set_label('normed FWHM',size=20,labelpad=10)
plt.xlim(0.1,9.3)
plt.ylim(0.1,9.3)
plt.ylabel('y-axis [x1000 pixels]',size=25,labelpad=15)
plt.xlabel('x-axis [x1000 pixels]',size=25,labelpad=15)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid()


#####SPLUS
mat,v1,v2 = R.CC_numberdensity_contour_rangefixed(fw,rad2,0.2,0,1,0)

base=N.arange(1.8,2.5,0.03)
plt.figure(2,figsize = (17,10),dpi=70, facecolor='w', edgecolor='k')
plt.clf()

dos = plt.axes([0.35,0.215,0.2,0.7])
a1,a2,a3 = plt.hist(fw2,base,facecolor='blue',alpha=0.15,orientation='horizontal',normed=1)
a1,a2,a3 = plt.hist(fw2,base,histtype='step',color='blue',alpha=0.5,orientation='horizontal',normed=1)
plt.ylim(2.4,3.);plt.grid()
plt.xticks(fontsize=20);plt.xlabel('frequence',size=25,labelpad=15)
plt.xlim(0,3)
# ;plt.legend(['$\#$4000'],loc='upper right',fontsize=12)

uno = plt.axes([0.1,.215,0.25,0.7])
plt.contour(v1,v2,mat,500,linewidths=2)
plt.ylabel('normed FWHM',size=25,labelpad=15)
plt.xlabel('radial distance [arcmin]',size=25,labelpad=15)
plt.xticks(fontsize=20);plt.yticks(fontsize=20)
plt.ylim(1.85,2.4);plt.grid();plt.xlim(1,59)

tres = plt.axes([0.625,0.215,0.3,0.7]);plt.scatter(x/1000., y/1000.,s=120,c=fw, marker=u'o', cmap=cm.gist_rainbow_r, alpha=0.2,vmin=1.5,vmax=2.6)
cb = plt.colorbar(pad=0.,format='%.2f');cb.set_label('normed FWHM',size=20,labelpad=10)
plt.xlim(0.1,9.3);plt.ylim(0.1,9.3);plt.ylabel('y-axis [x1000 pixels]',size=25,labelpad=15);plt.xlabel('x-axis [x1000 pixels]',size=25,labelpad=15)
plt.xticks(fontsize=20);plt.yticks(fontsize=20);plt.grid()



