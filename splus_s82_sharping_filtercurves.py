__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N

root_to_bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root_to_filters = root_to_bpz+'FILTER/'
final_root = root_to_filters+'SPLUS_cutted/'
if not os.path.exists(final_root):
    cmd8 = '/bin/mkdir %s '%(final_root)
    os.system(cmd8)

filters = U.get_str(root_to_filters+'laura_splus2.list',0)
nf = len(filters)

for ii in range(nf):
    xf,yf = U.get_data(root_to_filters+filters[ii],(0,1))
    yfr = N.where(yf<1.0e-4,0.0000,yf)
    new_file_name = final_root+filters[ii][:-4]+'.cut.res'
    U.put_data(new_file_name,(xf,yfr),'# w T')



