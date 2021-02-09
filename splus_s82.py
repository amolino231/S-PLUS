__author__ = 'albertomolino'

#/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N

root2cats = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Sept2017/'
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root2codes = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/'
lista_cats = 'spz.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)
#splus_columns = root2cats + 'splus11.columns'
splus_columns = root2cats + 'splus11_petro.columns'
#splus_columns = root2cats + 'splus11_auto.columns'

for ii in range(n_cats):
    single_cat = cats_names[ii]
    bpz_cat = single_cat[:-3]+'bpz'
    if not os.path.exists(bpz_cat):
          flux_cat = single_cat[:-3]+'flux_comparison'
          cmd2  = 'python %sbpz.py %s '%(root2bpz,single_cat)
          cmd2 += '-COLUMNS %s -OUTPUT %s '%(splus_columns,bpz_cat)
          cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.01 -INTERP 5 -ZMAX 1.0 '
          cmd2 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(flux_cat)
          cmd2 += '-SPECTRA GOSMOSeB11.list -PRIOR flat -ZMIN 0.00001 -DZ 0.0001 '
          os.system(cmd2)


    # Remove temporal files created by BPZ.
    cmd3 = '/bin/rm %s*.npy'%(root2cats)
    os.system(cmd3)
