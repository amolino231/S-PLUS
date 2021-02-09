__author__ = 'albertomolino'

"""
This version differs from the previous one since now
it takes into account the color (g-r) distribution of
stars in SDSS/S82.
"""

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import coeio as C
import useful as U

root_to_cat = '/Users/albertomolino/doctorado/articulos/SPLUS/'
root_to_cat += 'calibration/sample4simulations/data/20180329/NGSL_convolved/'
catalog = root_to_cat + 'sample_final_models_convolved_04042018.cat'
catalog_sdss_colors = root_to_cat + 'IveReal_GR_normed.cat'
models = U.get_str(root_to_cat+'ngsl_pick_wd.list',0)
new_catalog = catalog[:-3]+'extended_gr.cat'

# display information?
verbose = 0

# Reading the general catalogue.
datos_main = C.loaddata(catalog)
head_main = C.loadheader(catalog)

# gr synthetic models
#gr_models = datos_main[:,6]-datos_main[:,8] #SPLUS
gr_models = datos_main[:,16]-datos_main[:,17] #SDSS

# Reading the SDSS_color catalog.
sdss_gr,density_gr = U.get_data(catalog_sdss_colors,(0,1))
n_color_bins = len(sdss_gr)-1

# Reading catalog dimensions.
n_stars_ori = N.shape(datos_main)[0]
n_cols = N.shape(datos_main)[1]
n_filters = n_cols - 1#12 #20  #####################################

# Magnitude range we want to expand.
m_min = 12  ####################################
m_max = 22.0  ####################################
delta_m = 0.01  ####################################
base_m = N.arange(m_min,m_max+delta_m,delta_m)
n_m_ele = len(base_m)

# Number of stars per magnitude-bin.
n_stars_rbin = 500  ####################

# Scaling the color distribution
final_density_gr = N.ceil(density_gr * n_stars_rbin)
new_n_stars_rbin = sum(final_density_gr) # becuase it may include a few more models.

# Final number of stars to simulate
n_stars_final = n_m_ele * int(new_n_stars_rbin)

# Final matrix where to store new data.
final_mat = N.zeros((n_stars_final,n_cols),float)

# defining rSDSS as the reference band.
ref_column = 8 #10 in this new catalogue r is at 8 ###################

#number additional variables (ex.: id, ra, dec,...)
nadd = 1 #3


# starting the loop...
kkk = 0
#n_m_ele = 1
for iii in range(n_m_ele):
    mr_reference = base_m[iii]
    print 'Creating a stellar sample with magnitude r = ',mr_reference

    for sss in range(n_color_bins):
        # here we select the color interval and
        # the amount of stars in each color interval
        good_color  = N.greater_equal(gr_models,sdss_gr[sss])
        good_color *= N.less_equal(gr_models,sdss_gr[sss+1])
        new_dim_gr_models = len(gr_models[good_color])
        print 'new_dim_gr_models',new_dim_gr_models
        print 'stars with %.2f<g-r<%.2f: %i '%(sdss_gr[sss],sdss_gr[sss+1],final_density_gr[sss])

        # Selecting those models from main catalogue.
        datos_redu = datos_main[good_color,:]
        #print N.shape(datos_redu)
        #pausa = raw_input('paused')

        # Selecting 'n_new_stars' random models.
        new_stars = N.random.random_integers(0,new_dim_gr_models-1,int(final_density_gr[sss]))
        n_new_stars = int(final_density_gr[sss])

        for jjj in range(n_new_stars):
            mr_orig_sim = datos_redu[new_stars[jjj],ref_column]
            dm = (mr_reference-mr_orig_sim)
            final_mat[kkk,:] = datos_redu[new_stars[jjj],:]
            final_mat[kkk,nadd:nadd+n_filters] += dm
            final_mat[kkk,0] = datos_redu[new_stars[jjj],0]

            if verbose:
                print 'new_star: ',new_stars[jjj]
                print 'mr_reference',mr_reference
                print 'mr_orig_sim',mr_orig_sim
                print 'dm: ',dm
                print 'before'
                print datos_main[new_stars[jjj],:]
                print 'after'
                print final_mat[kkk,:]

            kkk += 1

# Saving new file
C.savedata(final_mat[0:kkk,:],new_catalog, dir="",header=head_main)



