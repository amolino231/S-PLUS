__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import bpz_tools as B
import matplotlib.pyplot as plt
import numpy as N

# Roots and Paths
bpz_path = '/Users/albertomolino/codigos/bpz-1.99.2/'
filter_path = bpz_path+'FILTER/'
sed_path = bpz_path+'SED/'
ab_path = bpz_path+'AB/'
sed_lib = 'greisel.list'
#sed_lib = 'tauAV3Z.list'
#sed_lib = 'GOSMOSeB11.list'
final_root_plots = sed_path + sed_lib[:-5]+'delta_colorz_ugz/'
if not os.path.exists(final_root_plots):
   cmd = '/bin/mkdir %s '%(final_root_plots)
   os.system(cmd)

new_file_out = final_root_plots + 'analysis.txt'
fileout = open(new_file_out,'w')
fileout.write('# Model d_gz_T1 d_ug_T1 d_gz_T2 d_ug_T2 d_tot_T1  d_tot_T2 \n')

# Filters
filter_x1 = 'SPLUS_gSDSS'
filter_x2 = 'SPLUS_zSDSS'
filter_y1 = 'SPLUS_uJAVA'
filter_y2 = 'SPLUS_gSDSS'

#Models
sed_models = U.get_str(sed_path+sed_lib,0)
n_models = len(sed_models)

# Empirical ranges.
z_bins = [0.03,0.04,0.06,0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.22]
nz = len(z_bins) #10

"""
#Old Numbers
gz_t1  = [1.3,1.40,1.40,1.5,1.5,1.6,1.6,1.7,1.7,1.8]
ug_t1  = [1.6,1.7,1.7,1.7,1.7,1.7,1.7,1.8,1.80,1.9]
gz_t2  = [1.4,1.5,1.5,1.5,1.6,1.6,1.7,1.8,1.9,1.9]
ug_t2  = [1.8,1.8,1.9,1.9,2.0,2.0,2.0,2.0,2.0,2.1]
"""

# New values
gz_t1  = [1.17,1.24,1.28,1.28,1.30,1.35,1.50,1.52,1.56,1.65,1.70]
ug_t1  = [1.85,1.85,1.85,1.85,1.85,1.85,1.85,1.85,1.85,1.85,1.85]
gz_t2  = [1.25,1.25,1.25,1.45,1.50,1.60,1.70,1.70,1.75,1.80,2.00]
ug_t2  = [2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25]

for ii in range(n_models):
    # Reading AB files for each template and filter.
    ab_filter_x1_model = ab_path+sed_models[ii][:-3]+filter_x1+'.AB'
    z_x1_model, mag_x1_model = U.get_data(ab_filter_x1_model,(0,1))
    ab_filter_x1 = B.flux2mag(mag_x1_model)

    ab_filter_x2_model = ab_path+sed_models[ii][:-3]+filter_x2+'.AB'
    z_x2_model, mag_x2_model = U.get_data(ab_filter_x2_model,(0,1))
    ab_filter_x2 = B.flux2mag(mag_x2_model)

    ab_filter_y1_model = ab_path+sed_models[ii][:-3]+filter_y1+'.AB'
    z_y1_model, mag_y1_model = U.get_data(ab_filter_y1_model,(0,1))
    ab_filter_y1 = B.flux2mag(mag_y1_model)

    ab_filter_y2_model = ab_path+sed_models[ii][:-3]+filter_y2+'.AB'
    z_y2_model, mag_y2_model = U.get_data(ab_filter_y2_model,(0,1))
    ab_filter_y2 = B.flux2mag(mag_y2_model)

    for zz in range(nz):
        pos_z  = N.argmin(abs(z_x1_model-z_bins[zz]))
        #Template Colour definition
        color_templ_xaxis = (ab_filter_x1[pos_z]-ab_filter_x2[pos_z])
        color_templ_yaxis = (ab_filter_y1[pos_z]-ab_filter_y2[pos_z])
        if zz<1:
           #First template
           delta_mx_t1 = abs(color_templ_xaxis-gz_t1[zz])
           delta_my_t1 = abs(color_templ_yaxis-ug_t1[zz])
           # Second template.
           delta_mx_t2 = abs(color_templ_xaxis-gz_t2[zz])
           delta_my_t2 = abs(color_templ_yaxis-ug_t2[zz])
        else:
           #First template
           delta_mx_t1 += abs(color_templ_xaxis-gz_t1[zz])
           delta_my_t1 += abs(color_templ_yaxis-ug_t1[zz])
           # Second template.
           delta_mx_t2 += abs(color_templ_xaxis-gz_t2[zz])
           delta_my_t2 += abs(color_templ_yaxis-ug_t2[zz])

    fileout.write('%s  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  \n'%(sed_models[ii],
    delta_mx_t1,delta_my_t1,delta_mx_t2,delta_my_t2,delta_mx_t1+delta_my_t1,delta_mx_t2+delta_my_t2))
fileout.close()

