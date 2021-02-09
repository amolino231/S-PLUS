__author__ = 'albertomolino'

import os,sys
import numpy as N
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
#import alhambra_photools as A
import matplotlib.pyplot as plt
from matplotlib import cm

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root_to_cats = root+'splus_cats_NGSL/'
final_path = root + 'data_quality/astrometry/'
lista_cats = root_to_cats+'photometry.list'
cats = U.get_str(lista_cats,0)
n_cats = len(cats)
# Reference Catalogue.
sdss_s82_spz_cat = '/Users/albertomolino/Postdoc/T80S_Pipeline/targets/SDSS_S82/S82Stars.cat' #stars!!!
#sdss_s82_spz_cat = root_to_cats + 'S82_SDSS.cat' #specz galaxies!!!

# It creates (if it doesn't exist already) a folder where to save the new data.
coordinates = root_to_cats + 'coo/'
if not os.path.exists(coordinates):
   cmd = '/bin/mkdir %s'%(coordinates)
   os.system(cmd)

# Creating (if necessary) new files containing IDs,RA,DEC from each S-PLUS/S82 tile.
for ii in range(n_cats):
    field = os.path.basename(cats[ii])[9:-15]
    new_radec_filename = coordinates+'STRIPE82_%s_radec.coo'%(field)
    if not os.path.exists(new_radec_filename):
           new_file = open(new_radec_filename,'w')
           print 'Reading file %i out of %i'%(ii+1,n_cats)
           ids,ra,dec,x,y = U.get_data(cats[ii],(0,1,2,3,4))
           good_xy = N.greater_equal(x,2500) * N.less_equal(x,7000)
           good_xy *= N.greater_equal(y,2500) * N.less_equal(y,7000)
           ids,ra,dec,x,y = U.multicompress(good_xy,(ids,ra,dec,x,y))
           n_gals = len(ids)
           for ss in range(n_gals):
               #linea = '%s  %i  %.7f  %.7f  \n'%(field,ids[ss],ra[ss],dec[ss])
               linea = '%s  %i  %.7f  %.7f  %.1f  %.1f  \n'%(field,ids[ss],ra[ss],dec[ss],x[ss],y[ss])
               new_file.write(linea)
           new_file.close()

# It performs a match between S-PLUS coo and SDSS coo to
for ggg in range(n_cats):
    radec_match_cat = coordinates + os.path.basename(cats[ggg][:-3] + 'radec.match2sdss.cat')
    if not os.path.exists(radec_match_cat):
       print 'Cross-matching file %i out of %i'%(ggg+1,n_cats)
       field = os.path.basename(cats[ggg])[9:-15]
       new_radec_filename = coordinates+'STRIPE82_%s_radec.coo'%(field)
       cmd_cross_match  = "java -jar /Users/albertomolino/codigos/Stilts/stilts.jar "
       cmd_cross_match += "tmatch2 ifmt1=ascii ifmt2=ascii in1=%s "%(new_radec_filename)
       cmd_cross_match += "in2=%s out=%s ofmt=ascii matcher=sky values1='$3 $4' "%(sdss_s82_spz_cat, radec_match_cat)
       cmd_cross_match += "values2='$1 $2' params=1 join=1and2 find=best progress=log"
       os.system(cmd_cross_match)

# Once all catalogues have been matched, it creates a master catalogue of coordinates.
# 1. Makes a file list with all catalogues.
matched_cat_list = coordinates + 'STRIPE82-radec.matched.coo.list'
if not os.path.exists(matched_cat_list):
   cmd = '/bin/ls %s*radec.match2sdss.cat > %s'%(coordinates,matched_cat_list)
   print cmd
   os.system(cmd)

# 2. It makes a master catalogue with all common detections.
if os.path.exists(matched_cat_list):
   master_radec_coo_files = coordinates + 'master.STRIPE82-radec.matched.coo.cat'
   if not os.path.exists(master_radec_coo_files):
      #Files to be read
      matched_cats = U.get_str(matched_cat_list,0)
      n_match_cats = len(matched_cats)
      # Creating the final master catalogue.
      master_file = open(master_radec_coo_files,'w')
      #master_file.write('# Field IdSP RASP DecSP RASD DecSD zs r i sep \n')
      master_file.write('# Field IdSP RASP DecSP X Y RASD DecSD zs r i sep \n')
      for hh in range(n_match_cats):
          print 'Reading file %i out of %i'%(hh+1,n_match_cats)
          temp_file = matched_cats[hh]
          temp = open(temp_file,'r')
          datos = temp.read()
          datos = datos.split('\n')
          temp.close()
          n_ele = N.shape(datos)[0] - 1
          for jj in range(n_ele):
              master_file.write(datos[jj+1]+' \n')
      master_file.close()


## It makes the plot...
if os.path.exists(master_radec_coo_files):
   #Reading data
   ra_sp,dec_sp,ra_sd,dec_sd,r_sd = U.get_data(master_radec_coo_files,(2,3,6,7,9))
   delta_ra  = (ra_sp - ra_sd) * 3600/0.55
   delta_dec = (dec_sp - dec_sd) * 3600/0.55
   good_sample = N.greater_equal(r_sd,14) * N.less_equal(r_sd,21.5)
   delta_ra_r,delta_dec_r,r_sd_r = U.multicompress(good_sample,(delta_ra,delta_dec,r_sd))

   # Starting the figure
   res = 50
   plt.figure(1,figsize = (12,10),dpi=70, facecolor='w', edgecolor='k')
   plt.clf()
   #plt.plot(U.mean_robust(delta_dec[::res])-0.15,U.mean_robust(delta_ra[::res]),'rs',ms=20)
   plt.scatter(delta_dec[::res]-0.15,delta_ra[::res],s=100,c=r_sd[::res], marker=u'o',cmap=cm.PuOr,alpha=0.25,vmin=14.0,vmax=21.5)
   #plt.plot(U.mean_robust(delta_dec[::res])-0.15,U.mean_robust(delta_ra[::res]),'rs',ms=20)
   #plt.legend(['Mean'],loc='upper left',fontsize=30)
   cb = plt.colorbar(pad=0.,format='%.1f')
   cb.set_label('R-band Magnitude',size=25,labelpad=10)
   plt.xlim(-1.99,1.99)
   plt.ylim(-1.99,1.99)
   plt.ylabel('$\delta$$(RA)$ [pix]',size=25,labelpad=5)
   plt.xlabel('$\delta$$(Dec)$ [pix]',size=25,labelpad=5)
   plt.xticks(fontsize=20)
   plt.yticks(fontsize=20)
   plt.grid()
   print '%.2f,%.2f'%(U.mean_robust(delta_dec[::res]),U.mean_robust(delta_ra[::res]))
   plt.savefig('/Users/albertomolino/Desktop/example1.png',dpi=100)

