__author__ = 'albertomolino'

"""
This code differs from 'splus_s82_SEDs_colorz_evol.py' in the following:
rather than using fixed iso-contours from SDSS data, it creates redshift-dependent
contours based on the S-PLUS/SDSS-S82 spec-z sample.
This fact makes possible to overplot the SED color-tracks as a function
of redshift and so understand what regions of the color-space are not
fully covered by the library of SED models.

"""

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import bpz_tools as B
import matplotlib.pyplot as plt

# Survey to be used...
splus = 1

#Paths
root_to_cats = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
root_to_cats += 'S82/Dec2017/splus_cats_NGSL/'

bpz_path = '/Users/albertomolino/codigos/bpz-1.99.2/'
filter_path = bpz_path+'FILTER/'
sed_path = bpz_path+'SED/'
ab_path = bpz_path+'AB/'
#sed_lib = 'TGOSMOSeB11new_recal.list'
#sed_lib = 'GOSMOSeB11new.list'
#sed_lib = 'COSMOSeB11new_recal.list'
sed_lib = 'eB11.list'
final_root_plots = sed_path + sed_lib[:-5]+'_colorz/'
if not os.path.exists(final_root_plots):
   cmd = '/bin/mkdir %s '%(final_root_plots)
   os.system(cmd)


# Redshift range.
z_min = 0.01
z_max = 0.50
delta_z = 0.02
z_range = N.arange(z_min,z_max+delta_z,delta_z)
res_z = len(z_range)

# Filters
filter_x1 = 'SPLUS_gSDSS'
filter_x2 = 'SPLUS_zSDSS'
filter_y1 = 'SPLUS_uJAVA'
filter_y2 = 'SPLUS_gSDSS'


#Models
sed_models = U.get_str(sed_path+sed_lib,0)
n_models = len(sed_models)

# Reading the S-PLUS/S82 data (spec-z cat)
#bpz_cat = root_to_cats+'bpz/'+'master.STRIPE82_Photometry.m21.bpz'
if splus:
    photo_cat = root_to_cats+'cat/'+'master.STRIPE82_Photometry.m21.cat'
    u,du,g,dg,r,dr,i,z,dz,zs = U.get_data(photo_cat,(15,16,60,61,78,79,96,114,115,125))
else:
    photo_cat = '/Users/albertomolino/Postdoc/T80S_Pipeline/targets/SDSS_S82/'
    photo_cat += 'catalogues/stripe82_spz_extcorr.cat'
    u,g,r,i,z,zs = U.get_data(photo_cat,(4,6,8,10,12,3))

clean_data = N.less_equal(abs(u-g),5.) * N.less_equal(abs(g-z),5.)
#clean_data *= N.less_equal(du,0.1) * N.less_equal(dg,0.07) * N.less_equal(dz,0.07)
u,g,r,i,z,zs = U.multicompress(clean_data,(u,g,r,i,z,zs))

# Defining colours
color_x = g-z
color_y = u-g

# Estimating dimensionality
ab_file_example = ab_path+sed_models[0][:-3]+filter_x1+'.AB'
z_ab = U.get_data(ab_file_example,0)
nz = len(z_ab)

# Creating matrix where to save the data.
ab_filter_x1 = N.zeros((nz,n_models),float)
ab_filter_x2 = N.zeros((nz,n_models),float)
ab_filter_y1 = N.zeros((nz,n_models),float)
ab_filter_y2 = N.zeros((nz,n_models),float)


for ii in range(n_models):
    # Reading AB files for each template and filter.
    ab_filter_x1_model = ab_path+sed_models[ii][:-3]+filter_x1+'.AB'
    z_x1_model, mag_x1_model = U.get_data(ab_filter_x1_model,(0,1))
    ab_filter_x1[:,ii] = B.flux2mag(mag_x1_model)

    ab_filter_x2_model = ab_path+sed_models[ii][:-3]+filter_x2+'.AB'
    z_x2_model, mag_x2_model = U.get_data(ab_filter_x2_model,(0,1))
    ab_filter_x2[:,ii] = B.flux2mag(mag_x2_model)

    ab_filter_y1_model = ab_path+sed_models[ii][:-3]+filter_y1+'.AB'
    z_y1_model, mag_y1_model = U.get_data(ab_filter_y1_model,(0,1))
    ab_filter_y1[:,ii] = B.flux2mag(mag_y1_model)

    ab_filter_y2_model = ab_path+sed_models[ii][:-3]+filter_y2+'.AB'
    z_y2_model, mag_y2_model = U.get_data(ab_filter_y2_model,(0,1))
    ab_filter_y2[:,ii] = B.flux2mag(mag_y2_model)

for zz in range(res_z-1):
    # Select galaxies with that redshift.
    if zz<1:
       good_z = N.less_equal(zs,z_range[zz+1])
    else:
       good_z  = N.greater_equal(zs,z_range[zz])
       good_z *= N.less_equal(zs,z_range[zz+1])

    mean_z = U.mean_robust(zs[good_z])
    pos_z  = N.argmin(abs(mean_z-z_ab))

    #Template Colour definition
    color_templ_xaxis = (ab_filter_x1[pos_z,:]-ab_filter_x2[pos_z,:])
    color_templ_yaxis = (ab_filter_y1[pos_z,:]-ab_filter_y2[pos_z,:])

    # Plots
    plt.figure(1,figsize = (12,10),dpi=75, facecolor='w', edgecolor='k')
    plt.clf()
    plt.plot(color_x[good_z],color_y[good_z],'bo',alpha=0.2,ms=10)
    if mean_z<0.2:
        plt.xlim(0.,2.5)
        plt.ylim(0.,2.5)
    else:
        plt.xlim(0.5,3.5)
        plt.ylim(0.,3.5)
    plt.grid()
    plt.plot(color_templ_xaxis,color_templ_yaxis,'-s',ms=4,lw=1,color='red',alpha=0.8)
    for ii in range(n_models):
        plt.annotate('%i'%(ii+1),(color_templ_xaxis[ii]-0.0,color_templ_yaxis[ii]+0.0),
                 color='red',size=20)
    plt.xlabel('%s - %s'%(filter_x1,filter_x2),size=30)
    plt.ylabel('%s - %s'%(filter_y1,filter_y2),size=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    if zz<1:
        plt.legend(['$z<$$%.3f$'%(z_range[zz+1])],loc='upper left',numpoints=1,fontsize=30)
    else:
        plt.legend(['$%.3f$$<z<$$%.3f$'%(z_range[zz],z_range[zz+1])],loc='upper left',numpoints=1,fontsize=30)

    if splus:
        plot_filename = '%s%s_vs_%s%s.z%.2f.splus.png'%(filter_x1,filter_x2,filter_y1,filter_y2,mean_z)
    else:
        plot_filename = '%s%s_vs_%s%s.z%.2f.sdss.png'%(filter_x1,filter_x2,filter_y1,filter_y2,mean_z)
    #plt.savefig(final_root_plots+plot_filename,dpi=100)
    #pausa = raw_input('paused')

