__author__ = 'albertomolino'

import os,sys
import numpy as N
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import splus_calib_tools as sct
import matplotlib.pyplot as plt
from matplotlib import cm

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root_to_cats = root+'splus_cats/'
root_to_ind  = root+'splus_indiv_cats/'
final_path = root + 'data_quality/'
lista_cats = root_to_cats+'photometry.list'
cats = U.get_str(lista_cats,0)
n_cats = len(cats)

# Creating the file where to save the data.
file_out_name = final_path+'psf.Rband_centeref.cat'
file_out = open(file_out_name,'w')
file_out.write('# x y fwhm \n')

# Starting the game...
for sss in range(n_cats):
    field = os.path.basename(cats[sss])[9:-15]
    master_cat = root_to_cats + 'STRIPE82-%s_Photometry.cat'%(field) ##
    print 'reading catalog %i out of %i '%(sss+1,n_cats)
    x,y,fwhm,mr = U.get_data(master_cat,(3,4,8,84))
    fwhm_master,stars = sct.get_seeing_from_data_pro(fwhm,mr)
    # Reading FWHM from Individual R images.
    ind_catalog = root_to_ind+'sex_STRIPE82-%s_R_swp.cat'%(field) ##
    fwhm_indiv = U.get_data(ind_catalog,6)
    x_stars,y_stars,fw_stars = U.multicompress(stars,(x,y,fwhm_indiv))
    center_pos_x = N.less(abs(x_stars-N.mean(x_stars)),1000.)
    center_pos_y = N.less(abs(y_stars-N.mean(y_stars)),1000.)
    fw_stars_center = fw_stars[center_pos_x * center_pos_y]
    for ii in range(len(x_stars)):
        linea = '%i  %i '%(x_stars[ii],y_stars[ii])
        mean_value_center = U.mean_robust(fw_stars_center)
        #valor = (fw_stars[ii]*fw_stars[ii])-(mean_value_center*mean_value_center)
        valor = (fw_stars[ii]-mean_value_center)
        linea += '%.2f \n'%(valor)
        file_out.write(linea)
file_out.close()

if os.path.exists(file_out_name):
    xx,yy,ff = U.get_data(file_out_name,(0,1,2))

# Starting the figure
res = 1
plt.figure(1,figsize = (12,10),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
plt.scatter(xx[::100]+30000.,yy[::100]+30000.,s=300,c=ff[::100],
            marker=u'o',cmap=cm.PuOr,alpha=0.9,vmin=-0.3,vmax=0.3)
cb = plt.colorbar(pad=0.,format='%.2f')
cb.set_label('FWHM - FWHM$_{center}$',size=25,labelpad=10)

plt.scatter(xx[::res]/1000.,yy[::res]/1000.,s=300,c=ff[::res],
            marker=u'o',cmap=cm.PuOr,alpha=0.15,vmin=-0.3,vmax=0.3)

#plt.legend(['Mean'],loc='upper left',fontsize=30)
#cb = plt.colorbar(pad=0.,format='%.2f')
#cb.set_label('FWHM - FWHM$_{center}$',size=25,labelpad=10)
plt.xlim(1,11)
plt.ylim(8.4,10.7)
plt.ylabel('Y-axis [10$^{3}$pix]',size=25,labelpad=1)
plt.xlabel('X-axis [10$^{3}$pix]',size=25,labelpad=2)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.grid()

#plt.savefig('/Users/albertomolino/Desktop/example1.png',dpi=100)



