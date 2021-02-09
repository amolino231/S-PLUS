__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as np
import numpy as U

root2models = '/Users/albertomolino/codigos/bpz-1.99.2/SED/NGSL/'
stellar_models = root2models+'NGSLnewsorted.list'

# Reading stellar models
models = U.get_str(stellar_models,0)
n_models = len(models)

# Intervals for new models
min_lambda = 1700.
max_lambda = 10000.
res_lambda = 2
interp_factor = 5

# New wavelength base for models
ll_base = np.arange(min_lambda,max_lambda+res_lambda,res_lambda)
dim_ll_base = len(ll_base)

#New dimension of models
new_mod_dim = n_models*interp_factor-1

# Here it creates a matrix with all models on the new ll_base
orig_mat = np.zeros((n_models,dim_ll_base),float)
new_mat  = np.zeros((new_mod_dim,dim_ll_base),float)

# Filling up the matrix
for ii in range(n_models):
    x,y = U.get_data(root2models+models[ii],(0,1))
    y_interp = np.interp(ll_base,x,y)
    orig_mat[ii,:] = y_interp[:]


# Interpolating models and storing data in new_mat matrix.
for ii in range(new_mod_dim):
    y1 = orig_mat[ii,:]
    y2 = orig_mat[ii+1,:]
    #for jj in range(interp_factor):



