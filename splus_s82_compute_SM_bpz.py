__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U

# General roots.
root = '/Volumes/CLASH/S82/specz/'
root_final = root + 'SM/'
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
# Files
cat_list = root + 'cats.list'
global_cats = U.get_str(cat_list,0)
n_cats = len(global_cats)

for ii in range(n_cats):
    cat  = global_cats[ii]
    cols = cat[:-3]+'cali.columns'
    bpz  = root_final + os.path.basename(cat)[:-3]+'bpz'
    if not os.path.exists(bpz):
       cmd  = 'python %sbpz.py %s '%(root2bpz,cat)
       cmd += '-COLUMNS %s -OUTPUT %s -CHECK yes '%(cols,bpz)
       cmd += '-SIGMA_EXPECTED 0.015 -INTERP 5 -FLUX_COMPARISON yes '
       cmd += '-SPECTRA COSMOSeB11new_recal.list -ZMIN 0.005 -PRIOR SM '
       cmd += '-DZ 0.001 -ZMAX 1.0 -HDF5 no -STELLAR_MASS yes -ONLY_TYPE yes '
       cmd += '-ABSOLUTE_MAGNITUDE yes ABSOLUTE_MAGNITUDE_FILTER SPLUS_rSDSS '
       cmd += '-ZP_ERRORS "0.03,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02" '
       print cmd
       os.system(cmd)

