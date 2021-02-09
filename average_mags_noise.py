__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import alhambra_photools as A
import matplotlib.pyplot as plt

"""
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Sept2017/'
root2cats = mainroot + 'data_quality/'
master_catalogue = root2cats + 'splus_master.cat'
master_columns   = root2cats + 'splus_auto.columns'
output_filename  = root2cats + 'magvsnoise.txt'
"""

master_catalogue = '/Users/albertomolino/jplus_data_download/SV02_March07/SV02_March07.clean.cat'
master_columns  = '/Users/albertomolino/jplus_data_download/SV02_March07/jplus11.columns'
output_filename = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Sept2017/data_quality/jplus_depth.txt'

# Reading data from master catalogue
mags  = A.get_magnitudes(master_catalogue,master_columns)
emags = A.get_errmagnitudes(master_catalogue,master_columns)

#defining variables
basem = N.arange(14,22.5,0.5)
values = N.zeros((len(basem),13),float)
values[:,0] = basem[:]

#starting the party
for ii in range(11):
    mm = mags[:,ii]
    em = emags[:,ii]
    good = N.greater_equal(mm,14) * N.less_equal(mm,22)
    mr = mm[good]
    emr = em[good]
    values[:,ii+1] = U.bin_stats(mr,emr,basem,'mean_robust')

U.put_2Darray(output_filename,values,'# mag emag_i','%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f ')




