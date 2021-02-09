__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U
import numpy as N
import bpz_tools as B
import phz_plots as P
import matplotlib.pyplot as plt

# FILES
root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
root += 'S82/Dec2017/splus_cats_NGSL/'
columns  = root + 'master_splus_auto.columns'
fluxcomp = root + 'COSMOSeB11new_recal_OT/masterOT_COSMOSeB11new_recal_SM.flux_comparison'

# starting.

filters = B.get_filter_list(columns)
nf = len(filters)

ft,fob,efob = P.getinfo4_mosaic_uncert_obsvsmod(columns,fluxcomp)
mo = U.get_data(fluxcomp,1)
g = N.greater_equal(mo,14) * N.less_equal(mo,20)

base = N.arange(14.,21.,1.)
basem_histo = N.arange(-1,1,0.01)
scatts = N.zeros(nf)
means  = N.zeros(nf)
for ii in range(nf):
    nameout = fluxcomp[:-15]+'%s.magzpdist.png'%(filters[ii][:-4])
    nameout2 = fluxcomp[:-15]+'magzpdist.txt'
    plt.figure(112, figsize = (8,10.),dpi=80, facecolor='w', edgecolor='k')
    plt.clf()
    mt = B.flux2mag(ft[ii])
    mob = B.flux2mag(fob[ii])
    emob = N.where(abs(fob[ii][:])>0,B.e_frac2mag(efob[ii][:])/fob[ii][:],-99)
    uno = plt.axes([0.15,.675,0.75,0.265])
    dm = mt[g]-mob[g]
    mmo = mo[g]
    emmob = emob[g]
    mtt = mt[g]
    mobb = mob[g]
    g2 = N.less(abs(dm),0.3)
    dm2,mo2,emo2,mt2,mob2 = U.multicompress(g2,(dm,mmo,emmob,mtt,mobb))
    a1,a2,a3=plt.hist(dm[g2],basem_histo,facecolor='grey',alpha=0.2)
    plt.xlim(-0.49,0.49)
    plt.grid()
    plt.title(filters[ii][:-4],size=20)
    # plt.xlabel('mt-mob',size=20)
    plt.ylabel('Counts',size=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    dos = plt.axes([0.15,.1,0.75,0.57])
    plt.plot(dm2[::10],mo2[::10],'bo',ms=2,alpha=0.1)
    plt.errorbar(dm2[::10],mo2[::10],[emo2[::10],emo2[::10]],fmt="bo",ms=2,alpha=0.1)
    linea = B.bin_stats(mo2,dm2,base,'mean_robust')
    plt.plot(linea,base,'r--',lw=5,alpha=0.5)
    plt.plot(linea*0.,base,'-',color='white',lw=3,alpha=0.9)
    plt.xlim(-0.49,0.49)
    plt.ylim(14,21)
    plt.xlabel('$\delta_{m}$',size=20)
    plt.ylabel('$r$',size=25)
    plt.grid()
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    means[ii]  = B.mean_robust(mt2-mob2)
    scatts[ii] = B.std_mad(mt2-mob2)
    print 'Filter %s: <offset> = %.3f, <scatter> = %.3f'%(filters[ii],means[ii],scatts[ii])
    plt.savefig(nameout,dpi=100)

    #deltas_ids_mo = nameout[:-3]+'ZPanals.txt'
    #U.put_data(deltas_ids_mo,(ids[g2]/1,mo[g2],zs[g2],tb[g2],dm2),'# IDS  MO  ZS  TB  Mt-Mob')
    # U.put_data(deltas_ids_mo,(dm[g2],ids[g2]/1,mo[g2]))
    # pausa = raw_input('paused')

# Saving dispersions
pepes = open(nameout2,'w')
pepes.write('#  Filter  Mean_Robust  Std_mad \n')
for ii in range(nf):
        ele = '%s  %.3f %.3f \n' %(filters[ii],means[ii], scatts[ii])
        pepes.write(ele)
pepes.close()
