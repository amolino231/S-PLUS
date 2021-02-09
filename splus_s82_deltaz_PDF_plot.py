__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import matplotlib.pyplot as plt
from splus_s82_photoz_precision_plot import master_PDZerrDistribution

#splus_s82_photoz_precision_plot.py
ruta  = '/Users/albertomolino/Postdoc/T80S_Pipeline/'
ruta += 'Commisioning/S82/Dec2017/splus_cats_NGSL/COSMOSeB11new_recal/'

hdf5list = ruta + 'PriorSM/HDF5/hdf5.list'
bpzlist  = ruta + 'PriorSM/bpzspz.list'
basez2b,basez3b,norm_dz_peaks_00,norm_dz_pdfs_00,ng_00 = master_PDZerrDistribution(hdf5list,bpzlist,21,0.0)
basez2b,basez3b,norm_dz_peaks_05,norm_dz_pdfs_05,ng_05 = master_PDZerrDistribution(hdf5list,bpzlist,21,0.5)
basez2b,basez3b,norm_dz_peaks_09,norm_dz_pdfs_09,ng_09 = master_PDZerrDistribution(hdf5list,bpzlist,21,0.9)

##
plt.figure(1,figsize = (13,6),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
plt.subplot(131)
plt.plot(basez3b,norm_dz_peaks_00/(norm_dz_peaks_00.max()),'-r',basez2b,norm_dz_pdfs_00/(1.*norm_dz_pdfs_00.max())-0.06,'b--',lw=6,alpha=0.75)
plt.xlim(-0.1,0.1)
plt.legend(['Peaks','PDFs'],loc='upper center',fontsize=22)
plt.grid()
plt.ylim(0.01,1.5)
# plt.title('$Odds$>0.0')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel('$\delta z$',size=24,labelpad=2)
# plt.ylabel('PDF',size=24,labelpad=10)
plt.ylabel('Probability  Density  Function',size=24,labelpad=10)
# ng_00 = 11729


plt.subplot(132)
plt.plot(basez3b,norm_dz_peaks_05/norm_dz_peaks_05.max(),'-r',basez2b,norm_dz_pdfs_05/(1.*norm_dz_pdfs_05.max())-0.06,'b--',lw=6,alpha=0.75)
plt.xlim(-0.1,0.1)
plt.ylim(0.01,1.5)
plt.legend(['Peaks','PDFs'],loc='upper center',fontsize=22)
plt.grid()
# plt.title('$Odds$>0.5')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel('$\delta z$',size=24,labelpad=2)
# ng_05 = 5290

plt.subplot(133)
plt.plot(basez3b,norm_dz_peaks_09/norm_dz_peaks_09.max(),'-r',basez2b,norm_dz_pdfs_09/norm_dz_pdfs_09.max(),'b--',lw=6,alpha=0.75)
plt.xlim(-0.1,0.1)
plt.legend(['Peaks','PDFs'],loc='upper center',fontsize=22)
plt.grid()
plt.ylim(0.01,1.5)
# plt.title('$Odds$>0.9')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel('$\delta z$',size=24,labelpad=2)
# ng_09 = 1071
