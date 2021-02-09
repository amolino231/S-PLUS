__author__ = 'albertomolino'

"""
The idea would be the following: combine the PDFs
from BPZ derived using 3 different apertures,
weighting each PDF by a normalized parameter:
Nodds_i = Oddi / Odd1+Odd2+Odd3

The question is the following:
1. would this method improve the reliability of PDFs?
2. would this method reduce the fraction of outliers? (bright,faint)

This idea has already be implemented in Carrasco-Kind+13 (TPZ) to
combine PDFs from different photo-z codes but not from different photometries!
They used a similar parameter to the Odds named "zConf" as a weight factor.

"""

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
#import bpz_tools as B

mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Sept2017/'
root2cats = mainroot + 'splus_cats/'
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root2codes = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/'
lista_cats = 'photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)
apertures = ['auto','petro','aper']
n_apers = len(apertures)
# splus_columns = root2cats + 'splus11_auto.columns'
sdss_s82_spz_cat = root2cats + 'S82_SDSS.cat'
# spectra = 'eB11.list'
spectra = 'GOSMOSeB11.list'

for ggg in range(n_cats):
    bpz_cali_1  = cats_names[ggg][:-3] + 'spz.%s_cali.bpz'%(apertures[0])
    bpz_cali_2  = cats_names[ggg][:-3] + 'spz.%s_cali.bpz'%(apertures[1])
    bpz_cali_3  = cats_names[ggg][:-3] + 'spz.%s_cali.bpz'%(apertures[2])