__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import numpy as np
import useful as U
import matplotlib.pyplot as plt

colores = ['blue','green','red','purple']

# Root to data
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2photo_depth = mainroot + 'data_quality/depth/'

bpz_cat = root2cats + 'COSMOSeB11new_recal/PriorSM/masterBPZ_PriorSM.bpz'


#root2files = '/Users/albertomolino/Postdoc/JPAS/PF/Data_CEFCA_Oct11/'
#final_root_cats = root2files + 'catalogues/reformatted/'
#root_to_complet = final_root_cats + 'completeness/'
#root_to_counts  = final_root_cats + 'counts/'
#final_root = final_root_cats + 'photoz_depth/'

''

# Files to be used.
aperture = 'MAG_RESTRICTED_PSFCOR'
complet_file_mr  = root_to_complet + 'AEGIS.master.spz.%s.'%(aperture)
complet_file_mr += 'cali.completeness.mr.cat'
number_counts = root_to_counts + 'AEGIS.master.z.counts.mr.cat'

# Reading the number counts.
base_mr,mr = U.get_data(number_counts,(0,1))
# Reading the completeness fraction
base_mo,cm_o1,cm_o2,cm_o3,cm_o4 = U.get_data(complet_file_mr,(0,1,2,3,4))

# Intervals to be used.
mmin = 17.
mmax = 28.
dm = 0.5
new_base = np.arange(mmin,mmax+dm,dm)

new_numb_counts = U.match_resol(base_mr,mr,new_base)
new_complet_o1 = U.match_resol(base_mo,cm_o1,new_base)
new_complet_o1 = np.where(new_complet_o1<0.,0.,new_complet_o1)
new_complet_o2 = U.match_resol(base_mo,cm_o2,new_base)
new_complet_o2 = np.where(new_complet_o2<0.,0.,new_complet_o2)
new_complet_o3 = U.match_resol(base_mo,cm_o3,new_base)
new_complet_o3 = np.where(new_complet_o3<0.,0.,new_complet_o3)
new_complet_o4 = U.match_resol(base_mo,cm_o4,new_base)
new_complet_o4 = np.where(new_complet_o4<0.,0.,new_complet_o4)

# Final number counts.
final_num_counts_o1 = new_numb_counts * new_complet_o1
final_num_counts_o2 = new_numb_counts * new_complet_o2
final_num_counts_o3 = new_numb_counts * new_complet_o3
final_num_counts_o4 = new_numb_counts * new_complet_o4

plt.figure(10,figsize=(10,10),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
plt.semilogy(new_base,final_num_counts_o1,color=colores[0],lw=8,alpha=0.5)
plt.semilogy(new_base,final_num_counts_o2,color=colores[1],lw=8,alpha=0.5)
plt.semilogy(new_base,final_num_counts_o3,color=colores[2],lw=8,alpha=0.5)
plt.semilogy(new_base,final_num_counts_o4,color=colores[3],lw=8,alpha=0.5)
plt.semilogy(new_base,new_numb_counts,'k--',lw=10,alpha=0.8)
plt.grid()
plt.legend(['$dz/(1+z)$$=$$0.010$',
            '$dz/(1+z)$$=$$0.007$',
            '$dz/(1+z)$$=$$0.005$',
            '$dz/(1+z)$$=$$0.003$',
            '$dz/(1+z)$$>$$0.010$'],
            loc='upper left',fontsize=22)
plt.xlim(17,23)
plt.ylim(10,1.0e+4)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.xlabel('$r$',size=30,labelpad=-1)
plt.ylabel('n($r$) [deg$^{-1}$]',size=30,labelpad=4)
plt.savefig(final_root+'AEGIS.master.spz.%s.cali.completeness.mr.1.png'%(aperture),dpi=100)

plt.figure(11,figsize=(10,10),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
plt.semilogy(new_base,np.cumsum(final_num_counts_o1),color=colores[0],lw=8,alpha=0.5)
plt.semilogy(new_base,np.cumsum(final_num_counts_o2),color=colores[1],lw=8,alpha=0.5)
plt.semilogy(new_base,np.cumsum(final_num_counts_o3),color=colores[2],lw=8,alpha=0.5)
plt.semilogy(new_base,np.cumsum(final_num_counts_o4),color=colores[3],lw=8,alpha=0.5)
plt.semilogy(new_base,np.cumsum(new_numb_counts),'k--',lw=10,alpha=0.8)
plt.grid()

plt.legend(['$dz/(1+z)$$=$$0.010$',
            '$dz/(1+z)$$=$$0.007$',
            '$dz/(1+z)$$=$$0.005$',
            '$dz/(1+z)$$=$$0.003$',
            '$dz/(1+z)$$>$$0.010$'],
            fontsize=22,loc='upper left')
plt.xlim(17,23)
plt.ylim(10,2.0e+4)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.xlabel('$r$',size=30,labelpad=-1)
plt.ylabel('Cumulative n($r$) [deg$^{-1}$]',size=30,labelpad=4)
plt.savefig(final_root+'AEGIS.master.spz.%s.cali.completeness.mr.2.png'%(aperture),dpi=100)


