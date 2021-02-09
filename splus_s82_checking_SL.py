__author__ = 'albertomolino'

#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

import os,sys
import numpy as N
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import splus_calib_tools as sct
#import matplotlib.pyplot as plt


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
mu_global = []
mg_global = []
mr_global = []
mi_global = []
mz_global = []
for ii in range(n_cats):
    catalog = cats[ii]
    print '%i out of %i'%(ii+1,n_cats)
    x,y,s2n,Flag,fw = U.get_data(catalog,(3,4,6,7,8))
    mu,mg,mr,mi,mz = U.get_data(catalog,(21,66,84,102,120)) #117,114
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
        mu_global.append(mu[kk])
        mg_global.append(mg[kk])
        mr_global.append(mr[kk])
        mi_global.append(mi[kk])
        mz_global.append(mz[kk])

n_total = len(x_global)
file_out_name = final_path+'ugriz_dividedbyseeing.cat'
file_out = open(file_out_name,'w')
file_out.write('# x y s2n Flag fwhm mu_aper mg_aper mr_aper mi_aper mz_aper \n')
for ss in range(n_total):
    linea = '%i  %i  %.2f  %i  %.3f   %.3f   %.3f  %.3f  %.3f  %.3f \n'%(x_global[ss],
                            y_global[ss],s2n_global[ss],Flag_global[ss],fw_global[ss],
                    mu_global[ss],mg_global[ss],mr_global[ss],mi_global[ss],mz_global[ss])
    file_out.write(linea)
file_out.close()
