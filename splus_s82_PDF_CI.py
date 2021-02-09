__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U
import numpy as N
import splus_s82_hdf5_tools as to
import matplotlib.pyplot as plt

# Roots and paths
sigma = '0150'
root = '/Volumes/CLASH/S82/specz/'

# Reading photometric catalogues.
bpz_cat_list = root + 'bpz_short.list'
bpz_cats = U.get_str(bpz_cat_list,0)
n_cats = len(bpz_cats)

check_mags  = 1
check_types = 1
check_odds  = 1

# Definition of variables.
CIS = []
MOS = []
ZSS = []
ODS = []
TBB = []
for ii in range(n_cats):
    hdf5_file = os.path.dirname(bpz_cats[ii])+'/HDF5/'
    hdf5_file += os.path.basename(bpz_cats[ii])[:-8]+'hdf5'

    if os.path.exists(hdf5_file):
       #print  hdf5_file
       print 'Reading file %i/%i '%(ii+1,n_cats)
       ids,zb,mo,odd,zs,tb = U.get_data(bpz_cats[ii],(0,9,10,5,9,4))
       n_gals = len(zs)
       #ci_value = to.compute_Wittman_CIs(hdf5_file,zs)
       gauss_factor = float('0.%s'%(sigma))
       ci_value = to.compute_Wittman_CIs_gaussConv(hdf5_file,zs,gauss_factor)

    for ss in range(n_gals):
        CIS.append(ci_value[ss])
        MOS.append(mo[ss])
        ZSS.append(zs[ss])
        ODS.append(odd[ss])
        TBB.append(tb[ss])


CIS2 = N.asarray(CIS)
MOS2 = N.asarray(MOS)
ZSS2 = N.asarray(ZSS)
ODS2 = N.asarray(ODS)
TBB2 = N.asarray(TBB)

baseci = N.arange(0.,1.01,0.01)
baseci2 = baseci[:-1]+((baseci[1]-baseci[0])/2.)
if check_mags:
    good1 = N.less_equal(MOS2,17)
    good2 = N.less_equal(MOS2,18)
    good3 = N.less_equal(MOS2,19)
    good4 = N.less_equal(MOS2,20)
    good5 = N.less_equal(MOS2,21)

    plt.figure(11)
    v1,v2,v3 = plt.hist(CIS2[good1],baseci,color='blue',alpha=0.3,cumulative=1)
    w1,w2,w3 = plt.hist(CIS2[good2],baseci,color='blue',alpha=0.3,cumulative=1)
    y1,y2,y3 = plt.hist(CIS2[good3],baseci,color='blue',alpha=0.3,cumulative=1)
    s1,s2,s3 = plt.hist(CIS2[good4],baseci,color='blue',alpha=0.3,cumulative=1)
    r1,r2,r3 = plt.hist(CIS2[good5],baseci,color='blue',alpha=0.3,cumulative=1)
    plt.figure(1)
    plt.clf()
    plt.plot(baseci2,v1/v1.max(),'r-')
    plt.plot(baseci2,w1/w1.max(),'b-')
    plt.plot(baseci2,y1/y1.max(),'g-')
    plt.plot(baseci2,s1/s1.max(),'c-')
    plt.plot(baseci2,r1/r1.max(),'m-')
    label1='$r<17:$ $%.1f$'%(abs(baseci2-(v1/v1.max())).sum())
    label2='$r<18:$ $%.1f$'%(abs(baseci2-(w1/w1.max())).sum())
    label3='$r<19:$ $%.1f$'%(abs(baseci2-(y1/y1.max())).sum())
    label4='$r<20:$ $%.1f$'%(abs(baseci2-(s1/s1.max())).sum())
    label5='$r<21:$ $%.1f$'%(abs(baseci2-(r1/r1.max())).sum())
    plt.legend([label1,label2,label3,label4,label5],loc='upper left',fontsize=25)
    final_plot_name = root+'SC/'+'Wittman.CGK0.%s.AB.png'%(sigma)
    plt.plot(baseci2,baseci2,'k--')
    plt.title('CONV_GAUSS_KERNEL = 0.%s'%(sigma),size=20)
    plt.grid()
    plt.xlabel('$C$',size=28,labelpad=3)
    plt.ylabel('$F(C)$',size=28,labelpad=3)
    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18)
    plt.savefig(final_plot_name,dpi=100)

if check_types:
    good1 = N.less_equal(TBB2,5.5) #tt[0:36] #red
    good2 = N.greater_equal(TBB2,5.5) #tt[36:] #blue

    plt.figure(11)
    v1,v2,v3 = plt.hist(CIS2[good1],baseci,color='red',alpha=0.3,cumulative=1)
    w1,w2,w3 = plt.hist(CIS2[good2],baseci,color='blue',alpha=0.3,cumulative=1)
    plt.figure(2)
    plt.clf()
    plt.plot(baseci2,v1/v1.max(),'r-')
    plt.plot(baseci2,w1/w1.max(),'b-')
    label1='$ET:$ $%.1f$'%(abs(baseci2-(v1/v1.max())).sum())
    label2='$LT:$ $%.1f$'%(abs(baseci2-(w1/w1.max())).sum())
    plt.legend([label1,label2],loc='upper left',fontsize=25)
    final_plot_name = root+'SC/'+'Wittman.CGK0.%s.SpT.png'%(sigma)
    plt.plot(baseci2,baseci2,'k--')
    plt.title('CONV_GAUSS_KERNEL = 0.%s'%(sigma),size=20)
    plt.grid()
    plt.xlabel('$C$',size=28,labelpad=3)
    plt.ylabel('$F(C)$',size=28,labelpad=3)
    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18)
    plt.savefig(final_plot_name,dpi=100)

if check_odds:
    good1 = N.greater_equal(ODS2,0.0)
    good2 = N.greater_equal(ODS2,0.3)
    good3 = N.greater_equal(ODS2,0.6)
    good4 = N.greater_equal(ODS2,0.9)

    plt.figure(11)
    v1,v2,v3 = plt.hist(CIS2[good1],baseci,color='blue',alpha=0.3,cumulative=1)
    w1,w2,w3 = plt.hist(CIS2[good2],baseci,color='blue',alpha=0.3,cumulative=1)
    y1,y2,y3 = plt.hist(CIS2[good3],baseci,color='blue',alpha=0.3,cumulative=1)
    s1,s2,s3 = plt.hist(CIS2[good4],baseci,color='blue',alpha=0.3,cumulative=1)
    plt.figure(3)
    plt.clf()
    plt.plot(baseci2,v1/v1.max(),'r-')
    plt.plot(baseci2,w1/w1.max(),'b-')
    plt.plot(baseci2,y1/y1.max(),'g-')
    plt.plot(baseci2,s1/s1.max(),'c-')
    label1='$O>0.0:$ $%.1f$'%(abs(baseci2-(v1/v1.max())).sum())
    label2='$O>0.3:$ $%.1f$'%(abs(baseci2-(w1/w1.max())).sum())
    label3='$O>0.6:$ $%.1f$'%(abs(baseci2-(y1/y1.max())).sum())
    label4='$O>0.9:$ $%.1f$'%(abs(baseci2-(s1/s1.max())).sum())
    plt.legend([label1,label2,label3,label4],loc='upper left',fontsize=25)
    final_plot_name = root+'SC/'+'Wittman.CGK0.%s.Odds.png'%(sigma)
    plt.plot(baseci2,baseci2,'k--')
    plt.title('CONV_GAUSS_KERNEL = 0.%s'%(sigma),size=20)
    plt.grid()
    plt.xlabel('$C$',size=28,labelpad=3)
    plt.ylabel('$F(C)$',size=28,labelpad=3)
    plt.yticks(fontsize=18)
    plt.xticks(fontsize=18)
    plt.savefig(final_plot_name,dpi=100)
