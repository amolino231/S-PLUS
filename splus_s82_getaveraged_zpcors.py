__author__ = 'albertomolino'

"""
This computes the averaged zero-point corrections
from a sample of columns files.

"""

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import alhambra_photools as A
import numpy as N

#Roots
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/COSMOSeB11new_recal/'
#bpz_cat = root2cats + 'master.STRIPE82_Photometry.m21_COSMOSeB11new_recal.bpz'

# List of cali.columns

list_cali_columns = root2cats + 'calicolumns.list'
if not os.path.exists(list_cali_columns):
   cmd = 'ls %s*.auto_cali.columns > %scalicolumns.list'%(root2cats,root2cats)
   os.system(cmd)


# Reading list
if os.path.exists(list_cali_columns):
   cols = U.get_str(list_cali_columns,0)
   n_c = len(cols)
   zp_val = N.zeros((12,n_c),float)
   zp_final = N.zeros(12)
   for ss in range(n_c):
       vars,evars,posref,zpe,zpc = A.get_usefulcolumns(cols[ss])
       #print zpc[:]
       #pausa = raw_input('paused')
       for ii in range(12):
           zp_val[ii,ss] = zpc[ii]

   for hh in range(12):
       zp_final[hh] = U.mean_robust(zp_val[hh,:])

   # Saving results
   master_cali_columns = root2cats + 'master_calicolumns.txt'
   U.put_data(master_cali_columns,(N.arange(12)+1,zp_final),'# Filter ZPc')

else:
    print 'File %s does not exist!'%(master_cali_columns)
