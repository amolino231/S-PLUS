__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import bpz_tools as B
import matplotlib.pyplot as plt

#Paths
bpz_path = '/Users/albertomolino/codigos/bpz-1.99.2/'
filter_path = bpz_path+'FILTER/'
sed_path = bpz_path+'SED/'
ab_path = bpz_path+'AB/'
sed_lib = 'GOSMOSeB11.list'
#sed_lib = 'GOSMOSeB11recal.list'
final_root_plots = sed_path + sed_lib[:-5]+'_colorz/'
if not os.path.exists(final_root_plots):
   cmd = '/bin/mkdir %s '%(final_root_plots)
   os.system(cmd)


# Redshift range.
z_min = 0.0
z_max = 0.2
delta_z = 0.01
z_range = N.arange(z_min,z_max+delta_z,delta_z)
res_z = len(z_range)

# Filters
filter_x1 = 'g_SDSS'
filter_x2 = 'z_SDSS'
filter_y1 = 'u_SDSS'
filter_y2 = 'g_SDSS'

#Models
sed_models = U.get_str(sed_path+sed_lib,0)
n_models = len(sed_models)

# SDSS/S82 Contours
sdss_s82_path = '/Users/albertomolino/Postdoc/T80S_Pipeline/targets/SDSS_S82/'
sdss_s82_cat = sdss_s82_path+'catalogues/stripe82_spz_extcorr.cat'

u_mag,g_mag,r_mag,i_mag,z_mag,z_spec = U.get_data(sdss_s82_cat,(4,6,8,10,12,3))
good_sample = N.less(abs(u_mag-g_mag),5.)
good_sample*= N.less(abs(g_mag-z_mag),5.)
good_sample*= N.greater_equal(z_spec,z_min)
good_sample*= N.less_equal(z_spec,z_max)
u_mag,g_mag,r_mag  = U.multicompress(good_sample,(u_mag,g_mag,r_mag))
i_mag,z_mag,z_spec = U.multicompress(good_sample,(i_mag,z_mag,z_spec))

#Reduce sample size!
sdss_res = 3
u_mag,g_mag,r_mag = u_mag[::sdss_res],g_mag[::sdss_res],r_mag[::sdss_res]
i_mag,z_mag,z_spec = i_mag[::sdss_res],z_mag[::sdss_res],z_spec[::sdss_res]

#Reading contours
path_contour_map = sdss_s82_path+'catalogues/contours/'
x_cmap,y_cmap = U.get_data(path_contour_map+'S82cont2.dat',(0,1))
z_cmap = U.get_2Darray(path_contour_map+'S82cont.dat')


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

for zz in range(res_z):
    pos_z = N.argmin(abs(z_range[zz]-z_ab))
    z_eff = z_range[zz]

    #Template Colour definition
    color_templ_xaxis = (ab_filter_x1[pos_z,:]-ab_filter_x2[pos_z,:])
    color_templ_yaxis = (ab_filter_y1[pos_z,:]-ab_filter_y2[pos_z,:])

    # Plots
    plt.figure(1,figsize = (12,10),dpi=75, facecolor='w', edgecolor='k')
    plt.clf()
    plt.contour(x_cmap,y_cmap,z_cmap,15,linewidths=2)
    #plt.xlim(-0.5,2.5)
    #plt.ylim(-0.5,2.5)
    plt.xlim(1.,2.2)
    plt.ylim(1.,2.)
    plt.grid()
    plt.plot(color_templ_xaxis,color_templ_yaxis,'-s',ms=4,lw=1,color='black',alpha=0.8)
    for ii in range(n_models):
        plt.annotate('%i'%(ii+1),(color_templ_xaxis[ii]-0.0,color_templ_yaxis[ii]+0.0),
                 color='black',size=20)
    plt.xlabel('%s - %s'%(filter_x1,filter_x2),size=30)
    plt.ylabel('%s - %s'%(filter_y1,filter_y2),size=30)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    plt.legend(['$z$ $=$ $%s$'%(z_eff)],loc='upper left',numpoints=1,fontsize=30)
    plot_filename = '%s%s_vs_%s%s.z%.2f.png'%(filter_x1,filter_x2,filter_y1,filter_y2,z_eff)
    plt.savefig(final_root_plots+plot_filename,dpi=100)
    #pausa = raw_input('paused')
