__author__ = 'albertomolino'

"""
This routine creates several plots showing the computed
n(z) of galaxies from the catalogues we have in the S82.
--
1. A comparison between spec-z and photo-z data.
--
This version differs from 'splus_s82_HDF5_ndz.py' because
it only uses galaxies with spec-z information.
===
Once the final SED+Prior is setup, i need to rerun BPZ on
the spec-z catalogues and create the corresponding HDF5 files.

"""

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U
import numpy as N
import splus_s82_hdf5_tools as to
import matplotlib.pyplot as plt

# General roots.
#root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
#root += 'S82/Dec2017/splus_cats_NGSL/COSMOSeB11new_recal/PriorSM/'
root = '/Volumes/CLASH/S82/specz/'

# Reading HDF5 files.
root_to_hdf5 = root + 'HDF5/'
final_root_hdf5 = root_to_hdf5 + 'PDFs/'
if not os.path.exists(final_root_hdf5):
   cmd = '/bin/mkdir %s '%(final_root_hdf5)
   os.system(cmd)
#hdf5list = root_to_hdf5 + 'hdf5_short.list'
#hdf5_files = U.get_str(hdf5list,0)
#n_hdf5 = len(hdf5_files)

# Reading photometric catalogues.
cat_list = root + 'bpz_short.list'
global_cats = U.get_str(cat_list,0)
n_cats = len(global_cats)

# Checking dimensionality!
#if n_cats != n_hdf5:
#   print 'Dimension missmatch!'
#   sys.exit()

odds_cut = 1

# Magnitude bins.
#mag_bins = [16,18,20]
mag_bins = [17,19,21]
n_mags = len(mag_bins)
zbs_1 = []
zbs_2 = []
zbs_3 = []
# Starting loop.
for ii in range(n_mags):
    print 'Processing detections with R<%i: '%(mag_bins[ii])
    if odds_cut:
       final_pdf_file = final_root_hdf5 + 'master_R%i_PDF.specz.odds09.txt'%(mag_bins[ii])
    else:
       final_pdf_file = final_root_hdf5 + 'master_R%i_PDF.specz.txt'%(mag_bins[ii])
    if not os.path.exists(final_pdf_file):
       kk = 0
       for ss in range(n_cats):
           #1. Select potential stars with R<Ri
           if odds_cut:
              ids,zb,mo,odd = U.get_data(global_cats[ss],(0,9,10,5))
              p_gal = N.ones(len(ids))
              p_gal = N.where(odd>0.89,1.,0.)
           else:
              ids,zb,mo = U.get_data(global_cats[ss],(0,9,10))
              p_gal = N.ones(len(ids))
           # Reading P(z) from HDF5 file.
           hdf5_file = os.path.dirname(global_cats[ss])+'/HDF5/'
           hdf5_file += os.path.basename(global_cats[ss])[:-8]+'hdf5'
           if os.path.exists(hdf5_file):
              print 'reading file %i/%i:  %s'%(ss+1,n_cats,os.path.basename(global_cats[ss])[:-8]+'hdf5')
              zz,pr,pb,pg,ng = to.getPDF_by_mag_and_weights(hdf5_file,mag_bins[ii],p_gal)
              #master_PDZerrDistribution
              # Storing P(z)
              if ss<1:
                 base_z = zz
                 final_pdf_red = pr
                 final_pdf_blue = pb
                 final_pdf_global = pg
                 n_gals = ng
              else:
                 final_pdf_red += pr
                 final_pdf_blue += pb
                 final_pdf_global += pg
                 n_gals += ng

              if odds_cut:
                 gm = N.less_equal(mo,mag_bins[ii])
                 gm *= N.greater(odd,0.89)
              else:
                  gm = N.less_equal(mo,mag_bins[ii])

              idr,zbr = U.multicompress(gm,(ids,zb))
              #print 'len(idr)',len(idr)
              for ff in range(len(idr)):
                  if ii < 1:
                     zbs_1.append(zbr[ff])
                  elif ii==1:
                     zbs_2.append(zbr[ff])
                  else:
                     zbs_3.append(zbr[ff])
                  kk+=1


       #print 'kk,ng',kk,n_gals
       #final_pdf_red *= n_gals
       #final_pdf_blue *= n_gals
       #final_pdf_global *= n_gals
       # Saving P(z|R<Ri)
       U.put_data(final_pdf_file,(zz,final_pdf_red,final_pdf_blue,final_pdf_global),'# z Pr Pb Pg ')

zbs_1 = N.asarray(zbs_1)
zbs_2 = N.asarray(zbs_2)
zbs_3 = N.asarray(zbs_3)


## PLOTS ##

#basez  = N.arange(0.005,0.8,0.01) #N.linspace(0.,0.8,200)
#basez2 = basez[:-1]+(basez[1]-basez[0])/2.

#root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
#root += 'S82/Dec2017/splus_cats_NGSL/COSMOSeB11new_recal/PriorSM/HDF5/PDFs/'

root = '/Volumes/CLASH/S82/specz/HDF5/PDFs/'
if odds_cut:
   z1,p1 = U.get_data(root+'master_R17_PDF.specz.odds09.txt',(0,3))
   z2,p2 = U.get_data(root+'master_R19_PDF.specz.odds09.txt',(0,3))
   z3,p3 = U.get_data(root+'master_R21_PDF.specz.odds09.txt',(0,3))
else:
   z1,p1 = U.get_data(root+'master_R17_PDF.specz.txt',(0,3))
   z2,p2 = U.get_data(root+'master_R19_PDF.specz.txt',(0,3))
   z3,p3 = U.get_data(root+'master_R21_PDF.specz.txt',(0,3))

if odds_cut:
   dz = 0.003
else:
   dz = 0.007

basez_1 = N.arange(0.01,max(zbs_1)+dz,dz)
basez2_1 = basez_1[:-1]+((basez_1[1]-basez_1[0])/2.)
v1,v2,v3 = plt.hist(zbs_1,basez_1,color='blue',alpha=0.3,normed=0)
p1r = U.match_resol(z1,p1,basez2_1)

basez_2 = N.arange(0.01,max(zbs_2)+dz,dz)
basez2_2 = basez_2[:-1]+((basez_2[1]-basez_2[0])/2.)
w1,w2,w3 = plt.hist(zbs_2,basez_2,color='blue',alpha=0.3,normed=0)
p2r = U.match_resol(z2,p2,basez2_2)

basez_3 = N.arange(0.01,max(zbs_3)+dz,dz)
basez2_3 = basez_3[:-1]+((basez_3[1]-basez_3[0])/2.)
q1,q2,q3 = plt.hist(zbs_3,basez_3,color='blue',alpha=0.3,normed=0)
p3r = U.match_resol(z3,p3,basez2_3)

plt.figure(21)
plt.clf()
plt.subplot(131)
plt.loglog(basez2_1,v1,'r-',lw=5,alpha=0.8)
plt.loglog(basez2_1,len(zbs_1)*(p1r/p1r.sum()),'b-',lw=5,alpha=0.8)
plt.xlabel('$z$',size=29,labelpad=-3)
plt.ylabel('$n(z)$',size=26,labelpad=1)
plt.legend(['spec-z','PDFs'],loc = 'lower left',fontsize=20)
plt.title('$r<%i$'%(mag_bins[0]),size=19)
plt.xlim(0.01,1.)
plt.ylim(1,2000)
plt.grid()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.subplot(132)
plt.loglog(basez2_2,w1,'r-',lw=5,alpha=0.8)
plt.loglog(basez2_2,len(zbs_2)*(p2r/p2r.sum()),'b-',lw=5,alpha=0.8)
plt.xlabel('$z$',size=29,labelpad=-3)
plt.legend(['spec-z','PDFs'],loc = 'lower left',fontsize=20)
plt.title('$r<%i$'%(mag_bins[1]),size=19)
plt.xlim(0.01,1.)
plt.ylim(1,2000)
plt.grid()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

plt.subplot(133)
plt.loglog(basez2_3,q1,'r-',lw=5,alpha=0.8)
plt.loglog(basez2_3,len(zbs_3)*(p3r/p3r.sum()),'b-',lw=5,alpha=0.8)
plt.xlabel('$z$',size=29,labelpad=-3)
plt.legend(['spec-z','PDFs'],loc = 'lower left',fontsize=20)
plt.title('$r<%i$'%(mag_bins[2]),size=19)
plt.xlim(0.01,1.)
plt.ylim(1,2000)
plt.grid()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)


# ===
plt.figure(3)
plt.clf()
plt.subplot(131)

plt.plot(basez2[::res3],v1r[::res3]/(v1r[::res3].sum()) - (p1r[::res3]/(p1r[::res3].sum())),'k-',lw=2,alpha=0.5)
plt.xlabel('$z$',size=24,labelpad=-3)
plt.ylabel('$n(z_{spz})$ - $n(z_{PDF})$',size=22,labelpad=1)
plt.title('$r<%i$'%(mag_bins[0]),size=19)
plt.xlim(0.,0.3)
#plt.ylim(0.,0.15)
plt.grid()

plt.subplot(132)
plt.plot(basez2[::res3],(w1r[::res3]/w1r[::res3].sum()) - (p2r[::res3]/p2r[::res3].sum()),'k-',lw=2,alpha=0.5)
plt.xlabel('$z$',size=24,labelpad=-3)
plt.title('$r<%i$'%(mag_bins[1]),size=19)
plt.xlim(0.,0.6)
#plt.ylim(0.,0.06)
plt.grid()

plt.subplot(133)
plt.plot(basez2[::res3],q1r[::res3]/(q1r[::res3].sum()) - (p3r[::res3]/(p3r[::res3].sum())),'k-',lw=2,alpha=0.5)
plt.xlabel('$z$',size=24,labelpad=-3)
plt.title('$r<%i$'%(mag_bins[2]),size=19)
plt.xlim(0.,1.0)
#plt.ylim(0.,0.042)
plt.grid()


