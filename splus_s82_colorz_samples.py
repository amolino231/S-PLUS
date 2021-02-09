__author__ = 'albertomolino'

"""
This routine compares in a couple of ways,
that the spectroscopic sample of galaxies
used to characterize the precision of S-PLUS
photo-z are representative of the main sample.

"""

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt

# S-PLUS filter names
filters = ['SPLUS_uJAVA','SPLUS_F0378W','SPLUS_F0395W',
           'SPLUS_F0410W','SPLUS_F0430W','SPLUS_gSDSS',
           'SPLUS_F0515W','SPLUS_rSDSS','SPLUS_F0660W',
           'SPLUS_iSDSS','SPLUS_F0861W','SPLUS_zSDSS']

n_filters = len(filters)

# BPZ galaxy models
sed_models = ['sedfit_restframe_z02_1091','sedfit_restframe_z01_2361',
              'sedfit_restframe_z00_1716','sedfit_restframe_z00_1613',
              'Ell3_A_0','Ell4_A_0','Ell5_A_0','Ell6_A_0','Ell7_A_0',
              'S0_A_0','Sa_A_0','Sa_A_1','Sb_A_0','Sb_A_1','Sbc_B10',
              'Scd_B10','SB1_B10','SB2_B10','SB3_B10','SB11_A_0_l']

n_seds = len(sed_models)

# Paths
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root_to_cats = mainroot + 'released_catalogues/'
root_to_plots = mainroot + 'data_quality/photoz/color_space/'
root_to_bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root_to_AB = root_to_bpz + 'AB/'
root_to_sed = root_to_bpz + 'SED/'

# Define filters for the color-z tracks (i.e., x2-x1)
filter_x1 = filters[7]
filter_x2 = filters[5]
filter_y1 = filters[11]
filter_y2 = filters[9]

# Define the new redshift range and resolution.
zmin = 0.
zmax = 0.5
new_dz = 0.01
base_z = N.arange(zmin,zmax+new_dz,new_dz)
n_z = len(base_z)

# Creating matrix where to save new colors.
mat_zcolors = N.zeros((n_z,n_seds),float)

for ii in range(n_seds):

    sed_ab_filter_x1 = sed_models[ii]+'.'+filter_x1+'.AB'
    sed_ab_filter_x1_z,sed_ab_filter_x1_f = U.get_data(sed_ab_filter_x1,(0,1))

    sed_ab_filter_x2 = sed_models[ii]+'.'+filter_x2+'.AB'
    sed_ab_filter_x2_z,sed_ab_filter_x2_f = U.get_data(sed_ab_filter_x2,(0,1))

    sed_ab_filter_y1 = sed_models[ii]+'.'+filter_y1+'.AB'
    sed_ab_filter_y1_z,sed_ab_filter_y1_f = U.get_data(sed_ab_filter_y1,(0,1))

    sed_ab_filter_y2 = sed_models[ii]+'.'+filter_y2+'.AB'
    sed_ab_filter_y2_z,sed_ab_filter_y2_f = U.get_data(sed_ab_filter_y2,(0,1))


# Reading photometric data from catalogue.
#catalog = root_to_cats + ''
#m_x1,s2n_x1 = U.get_data(catalog,())
#m_x2,s2n_x2 = U.get_data(catalog,())
#m_y1,s2n_y1 = U.get_data(catalog,())
#m_y2,s2n_y2 = U.get_data(catalog,())

#clean_sample = N.less(abs(m_x1),20) *