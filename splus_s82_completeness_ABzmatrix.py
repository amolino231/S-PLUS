__author__ = 'albertomolino'

"""
Here we use the galaxy sample with spectroscopic redshifts to
calculate the completeness of our observations.
--- September 2018 ---
"""

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as np
#import splus_calib_tools as ct

plots = 1

mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/'
lista_cats = 'cat/photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)
sdss_s82_spz_cat = root2cats + 'S82_SDSS_z1.cat'

########## Cross-matching with SDSS/S82.
for ggg in range(n_cats):
    splus_missing_speczcat = cats_names[ggg][:-3] + 'missing_spz.cat'
    splus_detected_speczcat = cats_names[ggg][:-3] + 'detected_spz.cat'

    if not os.path.exists(splus_missing_speczcat) or \
            not os.path.exists(splus_detected_speczcat):
       # I need to create a small catalog
       # containing galaxies within a common area
       # based on S-PLUS catalogs.

       ra_p,dec_p = U.get_data(cats_names[ggg],(1,2))
       x_p,y_p = U.get_data(cats_names[ggg],(3,4))
       ra_s,dec_s,z_s,r_s = U.get_data(sdss_s82_spz_cat,(0,1,2,3))

       # Let's remove edges.
       xmin = np.min(x_p)+300
       xmax = np.max(x_p)-300
       ymin = np.min(y_p)+300
       ymax = np.max(y_p)-300
       good_photo  = np.greater_equal(x_p,xmin)
       good_photo *= np.less_equal(x_p,xmax)
       good_photo *= np.greater_equal(x_p,ymin)
       good_photo *= np.less_equal(x_p,ymax)
       # Compress the photometric sample.
       ra_p,dec_p = U.multicompress(good_photo,(ra_p,dec_p))

       # Common area
       area_ra  = np.greater_equal(ra_s,min(ra_p))
       area_ra *= np.less_equal(ra_s,max(ra_p))
       area_dec  = np.greater_equal(dec_s,min(dec_p))
       area_dec *= np.less_equal(dec_s,max(dec_p))
       good_area = area_ra * area_dec
       #Compressing
       ra_s,dec_s,z_s,r_s = U.multicompress(good_area,(ra_s,dec_s,z_s,r_s))
       # Saving
       redu_spec = cats_names[ggg][:-3] + 'spez.areacommon.cat'
       U.put_data(redu_spec,(ra_s,dec_s,z_s,r_s),'# ra dec zs mr')

       # Here we find the missing galaxies
       cmd_cross_match  = "java -jar /Users/albertomolino/codigos/Stilts/stilts.jar "
       cmd_cross_match += "tmatch2 ifmt1=ascii ifmt2=ascii in1=%s "%(redu_spec)
       cmd_cross_match += "in2=%s out=%s ofmt=ascii matcher=sky values1='$1 $2' "%(cats_names[ggg],splus_missing_speczcat)
       cmd_cross_match += "values2='$2 $3' params=3 join=1not2 find=best progress=log"
       os.system(cmd_cross_match)

       # Here we find the detected galaxies
       cmd_cross_match  = "java -jar /Users/albertomolino/codigos/Stilts/stilts.jar "
       cmd_cross_match += "tmatch2 ifmt1=ascii ifmt2=ascii in1=%s "%(redu_spec)
       cmd_cross_match += "in2=%s out=%s ofmt=ascii matcher=sky values1='$1 $2' "%(cats_names[ggg],splus_detected_speczcat)
       cmd_cross_match += "values2='$2 $3' params=3 join=1and2 find=best progress=log"
       os.system(cmd_cross_match)


# Here we create master catalogs for both samples.
new_root = mainroot+'splus_cats_NGSL/cat/'
lista_missing = new_root+'SPLUS_missing_spz.cat.list'
lista_detected = new_root+'SPLUS_detected_spz.cat.list'
# Execute commands.
cmd0 = '/bin/ls %sSPLUS*Photometry.missing_spz.cat > %s'%(new_root,lista_missing)
os.system(cmd0)
cmd1 = '/bin/ls %sSPLUS*Photometry.detected_spz.cat > %s'%(new_root,lista_detected)
os.system(cmd1)

# It creates the final matrices to work with.
master_missing_catalog = new_root + 'master_missing_specz.cat'
if not os.path.exists(master_missing_catalog):
   import splus_calib_tools as ct
   ct.appendlistcatalog(lista_missing,master_missing_catalog)
master_deteced_catalog = new_root + 'master_detected_specz.cat'
if not os.path.exists(master_deteced_catalog):
   import splus_calib_tools as ct
   ct.appendlistcatalog(lista_detected,master_deteced_catalog)

## Here it creates the matrix
# Defining some parameters.
mmin=14
mmax=21
dm=0.5
zmin = 0.
zmax = 1.
dz = 0.05
base_m = np.arange(mmin,mmax+dm,dm)
nm = len(base_m)
base_z = np.arange(zmin,zmax+dz,dz)
nz = len(base_z)
ABz_matrix = np.zeros((nz-1,nm-1),float)

# Reading master catalogs.
zs_m,mr_m = U.get_data(master_missing_catalog,(2,3))
zs_d,mr_d = U.get_data(master_deteced_catalog,(2,3))

for ii in range(nz-1):
    # missing
    z_bin_m  = np.greater_equal(zs_m,base_z[ii])
    z_bin_m *= np.less_equal(zs_m,base_z[ii+1])
    zs_m_redu,mr_m_redu = U.multicompress(z_bin_m,(zs_m,mr_m))
    # detected
    z_bin_d  = np.greater_equal(zs_d,base_z[ii])
    z_bin_d *= np.less_equal(zs_d,base_z[ii+1])
    zs_d_redu,mr_d_redu = U.multicompress(z_bin_d,(zs_d,mr_d))
    print '==============================='
    for jj in range(nm-1):
        # missing
        r_bin_m  = np.greater_equal(mr_m_redu,base_m[jj])
        r_bin_m *= np.less_equal(mr_m_redu,base_m[jj+1])
        zs_m_redu2,mr_m_redu2 = U.multicompress(r_bin_m,(zs_m_redu,mr_m_redu))
        # detected
        r_bin_d  = np.greater_equal(mr_d_redu,base_m[jj])
        r_bin_d *= np.less_equal(mr_d_redu,base_m[jj+1])
        zs_d_redu2,mr_d_redu2 = U.multicompress(r_bin_d,(zs_d_redu,mr_d_redu))

        # Filling the matrix
        if len(zs_d_redu2)>0:
           complet = len(zs_m_redu2)/float(len(zs_d_redu2)+len(zs_m_redu2))
           ABz_matrix[ii,jj] = 100. - (complet*100.)
        else:
           #print 'empty bin'
           ABz_matrix[ii,jj] = -1.
        #print 'len(zs_d_redu)',len(zs_d_redu2)
        #print 'len(zs_m_redu)',len(zs_m_redu2)
        linea = '%.2f < m < %.2f & '%(base_m[jj],base_m[jj+1])
        linea += '%.2f < z < %.2f : '%(base_z[ii],base_z[ii+1])
        linea += '%.2f (D:%i, M:%i)'%(complet,len(zs_d_redu2),len(zs_m_redu2))
        #print linea

# Saving the matrix
#final_name_matrix = new_root + 'completeness_zAB.mat'
#U.put_2Darray(final_name_matrix,ABz_matrix)


if plots:
   import matplotlib.pyplot as plt
   from matplotlib.pyplot import cm

   base_m2 = base_m[:-1]+(base_m[1]-base_m[0])/2.
   base_z2 = base_z[:-1]+(base_z[1]-base_z[0])/2.

   plt.figure(1,figsize=(9,8),dpi=80, facecolor='w', edgecolor='k')
   plt.clf()
   for ii in range(nm-2):
       for jj in range(nz-2):
           if ABz_matrix[jj,ii]>0:
              plt.scatter(base_z2[jj],base_m2[ii],s=400,
                    c=ABz_matrix[jj,ii]/100.,marker=u's',
                    cmap=cm.PuOr,
                    alpha=0.85,vmin=0.0,vmax=1.0)
              print ABz_matrix[jj,ii]/100.

   #cmap=cm.PuOr,alpha=0.85, cm.PuOr
   #cm.jet
   cb = plt.colorbar(pad=0.,format='%.1f',ticks=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
   cb.set_label('Observed Completeness [%]',size=22,labelpad=10)
   plt.grid()
   plt.xlim(0.,1.0)
   plt.ylim(mmin,mmax-dm)
   #plt.legend(label4legend,loc='lower right',fontsize=20)
   plt.xticks(fontsize=20)
   plt.yticks(fontsize=20)
   plt.ylabel('$r$',size=30,labelpad=4)
   plt.xlabel('$z$',size=35,labelpad=-14)

