__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import coeio as C
import useful as U

root_to_cat = '/Users/albertomolino/doctorado/articulos/SPLUS/'
root_to_cat += 'calibration/sample4simulations/data/20180329/NGSL_convolved/'
#catalog = root_to_cat + 'SPLUS_SDSS_filters_4Laura.list_sample_stars_'
#catalog += 'SDSSDR10_aa_chi2_red_min.eq.0.03_20102017.cat'
catalog = root_to_cat + 'sample_NGSL_convolved_29032018.dat'
new_catalog = catalog[:-3]+'extended.cat'

# display information?
verbose = 0

# Reading the catalogue.
datos = C.loaddata(catalog)
head = C.loadheader(catalog)

# Reading catalog dimensions.
n_stars_ori = N.shape(datos)[0]
n_cols = N.shape(datos)[1]
n_filters = 12 #20  #####################################

# Magnitude range we want to expand.
m_min = 12  ####################################
m_max = 21  ####################################
delta_m = 0.01  ####################################
base_m = N.arange(m_min,m_max+delta_m,delta_m)
n_m_ele = len(base_m)

# Number of stars per magnitude-bin.
n_stars_final = 200 ####################################

# Final matrix where to store new data.
final_mat = N.zeros((n_m_ele*n_stars_final,n_cols),float)

# defining rSDSS as the reference band.
ref_column = 8 #10 in this new catalogue r is at 8 ###################

# Reading models from list
models = U.get_str(root_to_cat+'NGSLsorted.list',0)
n_models = len(models)
# Reading models from catalogue.
cat_models = U.get_str(catalog,0)


#number additional variables (ex.: id, ra, dec,...)
nadd = 1 #3

# starting the loop...
kkk = 0
for iii in range(n_m_ele):
    mr_reference = base_m[iii]
    print 'Creating a stellar sample with magnitude r = ',mr_reference
    # Randomly selects "n_stars_final" stars from original sample.
    new_stars = N.random.random_integers(0,n_stars_ori-1,n_stars_final)

    for jjj in range(n_stars_final):
        # Identify what position that model corresponds
        # to the NGSLsorted.list list.
        for sss in range(n_models):
            if cat_models[new_stars[jjj]] == models[sss]:
               modelo = sss+1
               if verbose:
                  print 'model_cat: ',cat_models[new_stars[jjj]]
                  print 'model_lib ', models[sss]
                  raw_input('paused')

        mr_orig_sim = datos[new_stars[jjj],ref_column]
        dm = (mr_reference-mr_orig_sim)
        final_mat[kkk,:] = datos[new_stars[jjj],:]
        final_mat[kkk,nadd:nadd+n_filters] += dm
        final_mat[kkk,0] = modelo

        if verbose:
           print 'new_star: ',new_stars[jjj]
           print 'mr_reference',mr_reference
           print 'mr_orig_sim',mr_orig_sim
           print 'dm: ',dm
           print 'before'
           print datos[new_stars[jjj],:]
           print 'after'
           print final_mat[kkk,:]

        kkk += 1

# Saving new file
C.savedata(final_mat,new_catalog, dir="",header=head)



