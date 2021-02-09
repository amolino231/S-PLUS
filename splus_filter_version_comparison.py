__author__ = 'albertomolino'

"""
This code plots the new and old filter curves of the SPLUS
and estimates the difference in the transmission efficience
between both sets (per filter) as the normalized difference
of the integrated signal within a filter.
"""


import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')

import useful as U
import pandas as pd
import bpz_tools as B
import numpy as N
import matplotlib.pyplot as plt

root = '/Users/albertomolino/codigos/bpz-1.99.2/FILTER/'
names = pd.read_table(root+'SPLUS_September2018/SPLUS_201809.list')
filtros_old = U.get_str(root+'SPLUS_July2017/SPLUS.list',0)
filtros_new = U.get_str(root+'SPLUS_September2018/SPLUS_201809.list',0)

base    = N.arange(3200,10000,10)
values  = N.zeros(12)
eff_wav = N.zeros(12)
for ii in range(12):
    eff_wav[ii] = B.effective_wavelength(filtros_new[ii])

for ii in range(12):
    x_o,y_o = U.get_data(root+filtros_old[ii],(0,1))
    x_n,y_n = U.get_data(root+filtros_new[ii],(0,1))
    y_o_r = U.match_resol(x_o,y_o,base)
    y_n_r = U.match_resol(x_n,y_n,base)
    values[ii] = N.sum(y_o_r-y_n_r)/N.sum(y_o_r)


plt.figure(20)
plt.subplot(211)
for ii in range(12):
    x,y = U.get_data(root+filtros_old[ii],(0,1))
    plt.fill_between(x*1.,y,0,alpha=0.4,color='grey',lw=1)
    x,y = U.get_data(root+filtros_new[ii],(0,1))
    plt.fill_between(x*1.,y,0,alpha=0.4,color='grey',lw=1)
plt.xlim(3000,10000)
plt.ylim(0.,0.8)
plt.ylabel('Throughput',size=28,labelpad=8)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title('S-PLUS Filter Curves: Old vs New',size=25)
plt.grid()

plt.subplot(212)
plt.plot(eff_wav,100.*values,'-s',color='grey',lw=6,ms=10)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('Difference [%]',size=30)
plt.xlabel('Wavelength [$\AA$]',size=28,labelpad=1)
plt.grid()

