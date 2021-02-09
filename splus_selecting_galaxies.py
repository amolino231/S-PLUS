__author__ = 'albertomolino'

import os,sys
#import numpy as N
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
#import coeio as C
import useful as U
import splus_calib_tools as sct
#import matplotlib.pyplot as plt
#from matplotlib import cm

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root_to_cats = root+'splus_final_cats/photo_cats/'
final_path = root + 'data_quality/'
lista_cats = root_to_cats+'photometrym23.list'
cats = U.get_str(lista_cats,0)
n_cats = len(cats)

for sss in range(n_cats):
    bpz_cat = cats[sss][:-3] + 'bpz'
    print 'bpz_cat',bpz_cat
    fwhm,mr = U.get_data(cats[sss],(8,78))
    galaxies = sct.get_galaxies_from_data(fwhm,mr)
    bright_galaxies = galaxies * N.less_equal(mr,21.0)
    new_bpz_cat = cats[sss][:-7]+'galsm21.bpz'
    if not os.path.exists(new_bpz_cat):
       sct.compress_bpz_catalogues(bpz_cat,bright_galaxies,new_bpz_cat)

