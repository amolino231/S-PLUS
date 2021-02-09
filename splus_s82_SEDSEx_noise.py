__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt
import bpz_tools as B
import alhambra_photools as A

"""
The idea of this routine is to test either
photometric noise and dispersion among models and data
are on the same scale as a function of the magnitude.
--
If these plots don't fall close to 1, it may mean that
photometric errors were either under/over-estimated.

"""

#Roots
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
#root2cats = mainroot + 'splus_cats_NGSL/COSMOSeB11new_recal/'
root2cats = mainroot + 'splus_cats_NGSL/COSMOSeB11new_recal_OT_after_SEDrecal/'
#bpz = root2cats + 'master.STRIPE82_Photometry.m21_COSMOSeB11new_recal_redu.bpz'
#fluxcomp = root2cats + 'master.STRIPE82_Photometry.m21_COSMOSeB11new_recal.flux_comparison'
fluxcomp = root2cats + 'master.STRIPE82_Photometry.m21_COSMOSeB11new_recal.OT.flux_comparison'
columns = root2cats + 'splus_auto_no_ids.columns'

# Final folder
final_root = root2cats+'plots/' # New root
if not os.path.exists(final_root):
    cmd8 = '/bin/mkdir %s '%(final_root)
    os.system(cmd8)

# Settings
m_min = 12
m_max = 21
delta_m = 1.



def SEDSEx_noise_comparison222(columns,fluxcomp,m_min,m_max,delta_m):
    """
fluxcomp = 'master.flux_comparison'
columns  = 'master.columns'
SEDSEx_noise_comparison(columns,fluxcomp,14,21,0.3)
    """

    basem = N.arange(m_min,m_max+delta_m,delta_m)

    ft,fob,efobs = A.get_fluxes(columns,fluxcomp)
    filters = B.get_filter_list(columns)
    nf = len(filters)

    for ss in range(nf):
        nband = ss
        mt = B.flux2mag(ft[:][nband])
        mo = B.flux2mag(fob[:][nband])
        emo = B.e_frac2mag(efobs[:][nband])/(fob[:][nband])
        g = N.less(abs(mo),28)
        mor,mtr,emor = U.multicompress(g,(mo,mt,emo))
        dm = mtr-mor
        g2 = N.less(abs(dm),1.)
        mor,mtr,emor = U.multicompress(g2,(mor,mtr,emor))

        linea_sdm = U.bin_stats(mor,dm,basem,'std_mad') #SED
        linea_sem = U.bin_stats(mor,emor,basem,'std_mad') #SEx
        linea_sem = N.where(linea_sem<1.0e-4,0.0,linea_sem)
        linea_dm  = U.bin_stats(mor,dm,basem,'mean_robust') #SED
        linea_em  = U.bin_stats(mor,emor,basem,'mean_robust') #SEx
        linea2    = linea_dm+linea_sdm  #SED
        linea1    = linea_sem+linea_em  #SEx
        """
        plt.figure(1, figsize = (11,9.5),dpi=80, facecolor='w', edgecolor='k')
        plt.clf()
        plt.subplot(211)
        plt.semilogy(basem,abs(linea_sdm),'r-',basem,linea_sem,'k-',lw=12)
        plt.legend(['SED','SEx'],loc='lower left',fontsize=30)
        plt.title('dispersions. Filter %s'%(filters[ss]),size=30)
        plt.grid()
        plt.xlim(m_min-delta_m,m_max+delta_m)
        plt.ylim(0.0001,1.)
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        # plt.xlabel('magnitude',size=22)
        plt.ylabel('$\delta m$',size=22)

        plt.subplot(212)
        plt.semilogy(basem,abs(linea_dm),'r-',basem,linea_em,'k-',lw=12)
        #plt.legend(['SED','SEx'],loc='lower left',fontsize=30)
        plt.title('mean',size=30)
        plt.grid()
        plt.xlim(m_min-delta_m,m_max+delta_m)
        plt.ylim(0.001,1.)
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        plt.xlabel('magnitude',size=22)
        plt.ylabel('$\delta m$',size=22)
        plt.savefig(final_root+'SEDdm1_Filter_%s.png'%(filters[ss]))
        """
        plt.figure(1, figsize = (11,9.5),dpi=80, facecolor='w', edgecolor='k')
        plt.clf()
        plt.semilogy(basem,abs(linea_sdm),'r-',basem,linea_sem,'k-',lw=12)
        plt.semilogy(basem,abs(linea_dm),'r--',basem,linea_em,'k--',lw=12)
        plt.title('Filter %s'%(filters[ss]),size=30)
        plt.xlim(m_min-delta_m,m_max+delta_m)
        plt.ylim(0.0001,1.)
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        plt.ylabel('$\delta m$',size=22)
        plt.legend(['$\sigma_{SED}$','$\sigma_{SEx}$','$\mu_{SED}$','$\mu_{SEx}$'],loc='lower left',fontsize=30)
        plt.grid()
        plt.xlabel('magnitude',size=22)
        plt.savefig(final_root+'SEDdm1_Filter_%s.png'%(filters[ss]))

        plt.figure(2, figsize = (11,9.5),dpi=80, facecolor='w', edgecolor='k')
        plt.clf()
        plt.semilogy(basem,abs(linea2),'r-',basem,linea1,'k-',lw=12)
        plt.legend(['SED','SEx'],loc='lower left',fontsize=30)
        plt.title('mean+disp. Filter %s'%(filters[ss]),size=30)
        plt.grid()
        plt.xlim(m_min-delta_m,m_max+delta_m)
        plt.ylim(0.01,1.)
        plt.xticks(fontsize=22)
        plt.yticks(fontsize=22)
        plt.xlabel('magnitude',size=22)
        plt.ylabel('$\delta m$',size=22)
        plt.savefig(final_root+'SEDdm2_Filter_%s.png'%(filters[ss]))


SEDSEx_noise_comparison222(columns,fluxcomp,m_min,m_max,delta_m)


