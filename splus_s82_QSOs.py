__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U

mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/'
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root2codes = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/'

spectra = 'QSOsorted'

"""
lista_cats = 'cat/individuals/photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)
lista_columns = root2cats + 'PriorSM/masterCOLUMNS_PriorSM.list'
cols_names = U.get_str(lista_columns,0)
n_cols = len(cols_names)
"""
splus_st_cat = mainroot + 'released_catalogues/S82_master_ps.cat'
splus_columns = root2cats + 'master_splus_auto.columns'

final_root = root2cats+'QSOs/' # New root
if not os.path.exists(final_root):
    cmd8 = '/bin/mkdir %s '%(final_root)
    os.system(cmd8)

nickname = os.path.basename(splus_st_cat)
bpz_out  = final_root+'QSOs.bpz'
flux_out = final_root+'QSOs.flux_comparison'
hdf5_out = final_root+'QSOs.hdf5'
if not os.path.exists(bpz_out):
   cmd2  = 'python %sbpz.py %s '%(root2bpz,splus_st_cat)
   cmd2 += '-COLUMNS %s -OUTPUT %s '%(splus_columns,bpz_out)
   cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.015 -INTERP 5 -ZMAX 10.0 '
   cmd2 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(flux_out)
   cmd2 += '-SPECTRA %s.list -ZMIN 0.005 -PRIOR flat -DZ 0.001 '%(spectra)
   cmd2 += '-HDF5 no -ONLY_TYPE no -USE_Z_S no '
   cmd2 += '-ZP_ERRORS "0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02" '
   os.system(cmd2)