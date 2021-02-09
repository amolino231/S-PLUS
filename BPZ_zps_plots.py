__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt

mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats/'
lista_cats = 'photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)

#apertures = ['auto','petro','aper']
#n_apers = len(apertures)
aperture = 'auto'

# filter names
filters = ['u','J0378','J0395','J0410','J0430','g',
          'J0515','r','J0660','i','J0861','z']

base_filtros = N.arange(12)+1

final_zpo = N.zeros((n_cats,12),float)
final_zpe = N.zeros((n_cats,12),float)

for ii in range(n_cats):
    cali_columns = cats_names[ii][:-4]+'.spz.z05.%s_cali.columns'%(aperture)
    final_zpe[ii,:],final_zpo[ii,:] = U.get_data(cali_columns,(3,4),12)

pepe = U.sum(final_zpe,axis=0)/(1.*n_cats)
pepo = U.sum(final_zpo,axis=0)/(1.*n_cats)

plt.figure(1, figsize=(14,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
plt.errorbar(base_filtros,pepo-pepo[7],pepe,fmt="-rs",alpha=0.5,ms=10,lw=7)
for ii in range(n_cats):
    plt.plot(base_filtros,final_zpo[ii,:]-final_zpo[ii,7],'k.',alpha=0.2,ms=10)
plt.errorbar(base_filtros,pepo-pepo[7],pepe,fmt="-rs",alpha=0.5,ms=10,lw=7)
plt.grid()
plt.legend(['$<zp_{c}>$'],loc='upper right',fontsize=35)
plt.xticks(base_filtros,filters,fontsize=20)
plt.yticks(fontsize=20)
plt.xlim(0.5,12.5)
plt.ylim(-0.2,0.2)
plt.ylabel('$zp_{c}$ based on SED-fitting',size=25,labelpad=5)
plt.xlabel('filters',size=25)

plt.savefig('/Users/albertomolino/Desktop/zps_offsets.png',dpi=100)


####
filters = ['u','J0378','J0395','J0410','J0430','g',
          'J0515','r','J0660','i','J0861','z']

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
#root += 'splus_cats_NGSL/JPLUS_specz_sample_March2018/NoPrior/specz/'
#lista1 = root + 'eB11/Master_COLUMNS_NoPrior_JPLUS_eB11.list'
#lista2 = root + 'COSMOSeB11new_recal/Master_COLUMNS_NoPrior_JPLUS_COSMOSeB11new_recal.list'
#splus_cats_NGSL/COSMOSeB11new_recal/NoPrior/masterCOLUMNS.list
lista1 = root + 'splus_cats_NGSL/COSMOSeB11new_recal/NoPrior/masterCOLUMNS.list'
lista2 = root + 'splus_cats_NGSL/COSMOSeB11new_recal/PriorSM/masterCOLUMNS_PriorSM.list'

base_filtros = N.arange(12)+1
lista = ''
cats_names1 = U.get_str(lista1,0)
cats_names2 = U.get_str(lista2,0)
n_cats_1 = len(cats_names1)
n_cats_2 = len(cats_names2)

final_zpo_1 = N.zeros((n_cats_1,12),float)
final_zpe_1 = N.zeros((n_cats_1,12),float)
final_zpo_2 = N.zeros((n_cats_2,12),float)
final_zpe_2 = N.zeros((n_cats_2,12),float)

for ii in range(n_cats_1):
    final_zpe_1[ii,:],final_zpo_1[ii,:] = U.get_data(cats_names1[ii],(3,4),12)
for ii in range(n_cats_2):
    final_zpe_2[ii,:],final_zpo_2[ii,:] = U.get_data(cats_names2[ii],(3,4),12)

av_zp_off_1 = U.sum(final_zpo_1,axis=0)/(1.*n_cats_1)
av_zp_off_2 = U.sum(final_zpo_2,axis=0)/(1.*n_cats_2)

#av_zp_err_1 = U.sum(final_zpe_1,axis=0)/(1.*n_cats_1)
#av_zp_err_2 = U.sum(final_zpe_2,axis=0)/(1.*n_cats_2)
av_zp_err_1 = N.zeros(12)
av_zp_err_2 = N.zeros(12)

for ii in range(12):
    av_zp_err_1[ii] = U.std_mad(final_zpo_1[ii,:])
    av_zp_err_2[ii] = U.std_mad(final_zpo_2[ii,:])

plt.figure(2, figsize=(14,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
aaa = av_zp_off_1-av_zp_off_1[7]
bbb = av_zp_off_2-av_zp_off_2[7]
plt.errorbar(base_filtros,(av_zp_off_1-N.mean(av_zp_off_1))*100.,av_zp_err_1*100.,fmt="-s",color='blue',alpha=0.5,ms=10,lw=7)
plt.errorbar(base_filtros,(av_zp_off_2-N.mean(av_zp_off_2))*100.,av_zp_err_2*100.,fmt="-s",color='red',alpha=0.5,ms=10,lw=7)
#plt.errorbar(base_filtros,(av_zp_off_3-N.mean(av_zp_off_3))*100.,av_zp_err_3*100.,fmt="-s",color='green',alpha=0.5,ms=10,lw=7)

#plt.errorbar(base_filtros,av_zp_off_2-av_zp_off_2[7],av_zp_err_2,fmt="-bs",alpha=0.5,ms=10,lw=7)
#plt.plot(base_filtros,bbb-aaa,'-ko',alpha=0.5,ms=10,lw=7)
plt.grid()
#plt.legend(['eB11','COSMOS'],loc='upper right',fontsize=35)
#plt.legend(['FLAT','SM'],loc='upper right',fontsize=35)
plt.legend(['$<zp_{c}>$'],loc='best',fontsize=35)
#plt.legend(['JPLUS (COSMOS-eB11)'],loc='best',fontsize=35)
#plt.xticks(base_filtros,filters,fontsize=20)
plt.yticks(fontsize=22)
plt.xticks(N.arange(12)+1,filters,size=22)
plt.xlim(0.5,12.5)
plt.ylim(-19.,19.)
plt.ylabel('$zp_{c}$ based on SED-fitting [x$10^{-2}$]',size=24,labelpad=0)
plt.xlabel('Filter',size=23)
