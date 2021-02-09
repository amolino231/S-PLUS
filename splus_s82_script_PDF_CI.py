__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U
import numpy as N
import splus_s82_hdf5_tools as to
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

# Roots and paths
root = '/Volumes/CLASH/S82/specz/'

# Reading photometric catalogues.
bpz_cat_list = root + 'bpz_short.list'
bpz_cats = U.get_str(bpz_cat_list,0)
n_cats = len(bpz_cats)

# Plots
check_mags  = 0
check_types = 0
check_odds  = 0
check_reds  = 1

# range
sigmas = N.linspace(0.005,0.05,15)+0.001
#sigmas = N.linspace(0.005,0.05,15)+0.001
#sigmas = sigmas[1:]
n_sigma = len(sigmas)

# Definition of variables.
for hhh in range(n_sigma):
    CIS = []
    MOS = []
    ZSS = []
    ODS = []
    TBB = []
    for ii in range(n_cats):
        hdf5_file = os.path.dirname(bpz_cats[ii])+'/HDF5/'
        hdf5_file += os.path.basename(bpz_cats[ii])[:-8]+'hdf5'

        if os.path.exists(hdf5_file):
           print 'Reading file %i/%i '%(ii+1,n_cats)
           mo,odd,zs,tb = U.get_data(bpz_cats[ii],(10,5,9,4))
           n_gals = len(zs)
           gauss_factor = float('%.3f'%(sigmas[hhh]))
           #ci_value = to.compute_Wittman_CIs(hdf5_file,zs)
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
        plt.plot(baseci2,baseci2,'k--')
        plt.title('CONV_GAUSS_KERNEL = %s'%(gauss_factor),size=20)
        plt.grid()
        plt.xlabel('$C$',size=28,labelpad=3)
        plt.ylabel('$F(C)$',size=28,labelpad=3)
        plt.yticks(fontsize=18)
        plt.xticks(fontsize=18)
        final_plot_name = root+'SC/'+'Wittman.CGK%s.AB.png'%(gauss_factor)
        plt.savefig(final_plot_name,dpi=100)
        final_file_name = root+'SC/'+'Wittman.CGK%s.AB.txt'%(gauss_factor)
        U.put_data(final_file_name,(baseci2,v1,w1,y1,s1,r1),'# C F17 F18 F19 F20 F21')

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
        plt.plot(baseci2,baseci2,'k--')
        plt.title('CONV_GAUSS_KERNEL = %s'%(gauss_factor),size=20)
        plt.grid()
        plt.xlabel('$C$',size=28,labelpad=3)
        plt.ylabel('$F(C)$',size=28,labelpad=3)
        plt.yticks(fontsize=18)
        plt.xticks(fontsize=18)
        final_plot_name = root+'SC/'+'Wittman.CGK%s.SpT.png'%(gauss_factor)
        plt.savefig(final_plot_name,dpi=100)
        final_file_name = root+'SC/'+'Wittman.CGK%s.SpT.txt'%(gauss_factor)
        U.put_data(final_file_name,(baseci2,v1,w1),'# C FET FLT')

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
        plt.plot(baseci2,baseci2,'k--')
        plt.title('CONV_GAUSS_KERNEL = %s'%(gauss_factor),size=20)
        plt.grid()
        plt.xlabel('$C$',size=28,labelpad=3)
        plt.ylabel('$F(C)$',size=28,labelpad=3)
        plt.yticks(fontsize=18)
        plt.xticks(fontsize=18)
        final_plot_name = root+'SC/'+'Wittman.CGK%s.Odds.png'%(gauss_factor)
        plt.savefig(final_plot_name,dpi=100)
        final_file_name = root+'SC/'+'Wittman.CGK%s.Odds.txt'%(gauss_factor)
        U.put_data(final_file_name,(baseci2,v1,w1,y1,s1),'# C FO00 FO03 FO06 FO09')

    if check_reds:
        good1 = N.less_equal(ZSS2,0.2)
        good2 = N.less_equal(ZSS2,0.4)
        good3 = N.less_equal(ZSS2,0.6)

        plt.figure(11)
        v1,v2,v3 = plt.hist(CIS2[good1],baseci,color='blue',alpha=0.3,cumulative=1)
        w1,w2,w3 = plt.hist(CIS2[good2],baseci,color='blue',alpha=0.3,cumulative=1)
        y1,y2,y3 = plt.hist(CIS2[good3],baseci,color='blue',alpha=0.3,cumulative=1)
        plt.figure(3)
        plt.clf()
        plt.plot(baseci2,v1/v1.max(),'r-')
        plt.plot(baseci2,w1/w1.max(),'b-')
        plt.plot(baseci2,y1/y1.max(),'g-')
        label1='$z<0.2:$ $%.1f$'%(abs(baseci2-(v1/v1.max())).sum())
        label2='$z<0.4:$ $%.1f$'%(abs(baseci2-(w1/w1.max())).sum())
        label3='$z<0.6:$ $%.1f$'%(abs(baseci2-(y1/y1.max())).sum())
        plt.legend([label1,label2,label3],loc='upper left',fontsize=25)
        plt.plot(baseci2,baseci2,'k--')
        plt.title('CONV_GAUSS_KERNEL = %s'%(gauss_factor),size=20)
        plt.grid()
        plt.xlabel('$C$',size=28,labelpad=3)
        plt.ylabel('$F(C)$',size=28,labelpad=3)
        plt.yticks(fontsize=18)
        plt.xticks(fontsize=18)
        final_plot_name = root+'SC/'+'Wittman.CGK%s.zs.png'%(gauss_factor)
        plt.savefig(final_plot_name,dpi=100)
        final_file_name = root+'SC/'+'Wittman.CGK%s.zs.txt'%(gauss_factor)
        U.put_data(final_file_name,(baseci2,v1,w1,y1),'# C FZ02 FZ04 FZ06')



###### FIGURE FOR PAPER #####
splus_colors = list(cm.Spectral_r(N.linspace(0, 1, 12)))

# AB
plt.figure(1)
plt.clf()
final_file_name_0 = root+'SC/AB/files/Wittman.CGK0.000.AB.txt'
final_file_name_1 = root+'SC/AB/files/Wittman.CGK0.019.AB.txt'
cc,v0 = U.get_data(final_file_name_0,(0,4))
cc,v1 = U.get_data(final_file_name_1,(0,4))
plt.plot(cc,v0/v0.max(),'-',lw=7,color=splus_colors[0])
plt.plot(cc,v1/v1.max(),'--',lw=7,color=splus_colors[0])
plt.plot(cc,cc,'-',color='black',lw=2,alpha=0.5)
label1='$Original$'
label2='$Smoothed$'
#plt.legend([label1,label2],loc='upper left',fontsize=25)
plt.grid()
plt.xlabel('$C$',size=30,labelpad=3)
plt.ylabel('$\hat F(C)$',size=30,labelpad=3)
plt.yticks(fontsize=21)
plt.xticks(fontsize=21)
plt.savefig(final_file_name_1[:-3]+'png',dpi=100)

# SpT
plt.figure(2)
plt.clf()
final_file_name_0 = root+'SC/SpT/files/Wittman.CGK0.000.SpT.txt'
final_file_name_1 = root+'SC/SpT/files/Wittman.CGK0.019.SpT.txt'
cc,v0e,v0l = U.get_data(final_file_name_0,(0,1,2))
cc,v1e,v1l = U.get_data(final_file_name_1,(0,1,2))
plt.plot(cc,v0e/v0e.max(),'r-',cc,v0l/v0l.max(),'b-',lw=7,alpha=0.7)
plt.plot(cc,v1e/v1e.max(),'r--',cc,v1l/v1l.max(),'b--',lw=7,alpha=0.7)
plt.plot(cc,cc,'-',color='black',lw=2,alpha=0.5)
label1='$Original$'
label2='$Smoothed$'
#plt.legend([label1,label2],loc='upper left',fontsize=25)
plt.grid()
plt.xlabel('$C$',size=30,labelpad=3)
plt.ylabel('$\hat F(C)$',size=30,labelpad=3)
plt.yticks(fontsize=21)
plt.xticks(fontsize=21)
plt.savefig(final_file_name_1[:-3]+'png',dpi=100)

# Odds
plt.figure(3)
plt.clf()
final_file_name_0 = root+'SC/Odds/files/Wittman.CGK0.000.Odds.txt'
final_file_name_1 = root+'SC/Odds/files/Wittman.CGK0.019.Odds.txt'
cc,v0 = U.get_data(final_file_name_0,(0,2))
cc,v1 = U.get_data(final_file_name_1,(0,2))
plt.plot(cc,v0/v0.max(),'-',lw=7,color=splus_colors[0])
plt.plot(cc,v1/v1.max(),'--',lw=7,color=splus_colors[0])
plt.plot(cc,cc,'-',color='black',lw=2,alpha=0.5)
label1='$Original$'
label2='$Smoothed$'
#plt.legend([label1,label2],loc='upper left',fontsize=25)
plt.grid()
plt.xlabel('$C$',size=30,labelpad=3)
plt.ylabel('$\hat F(C)$',size=30,labelpad=3)
plt.yticks(fontsize=21)
plt.xticks(fontsize=21)
plt.savefig(final_file_name_1[:-3]+'png',dpi=100)

# zs
plt.figure(4)
plt.clf()
final_file_name_0 = root+'SC/zs/files/Wittman.CGK0.000.zs.txt'
final_file_name_1 = root+'SC/zs/files/Wittman.CGK0.019.zs.txt'
cc,v0 = U.get_data(final_file_name_0,(0,3))
cc,v1 = U.get_data(final_file_name_1,(0,3))
plt.plot(cc,v0/v0.max(),'-',lw=7,color=splus_colors[0])
plt.plot(cc,v1/v1.max(),'--',lw=7,color=splus_colors[0])
plt.plot(cc,cc,'-',color='black',lw=2,alpha=0.5)
label1='$Original$'
label2='$Smoothed$'
#plt.legend([label1,label2],loc='upper left',fontsize=25)
plt.grid()
plt.xlabel('$C$',size=30,labelpad=3)
plt.ylabel('$\hat F(C)$',size=30,labelpad=3)
plt.yticks(fontsize=21)
plt.xticks(fontsize=21)
plt.savefig(final_file_name_1[:-3]+'png',dpi=100)



