__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as np
import matplotlib.pyplot as plt

root  = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
root += 'S82/Dec2017/data_quality/depth/'

cat = root + 'mr_petro_gals.cat'
outfilename = cat[:-3]+'depth.cat'

total_area = 1./80.

mmin = 14.
mmax = 22.
dm = 0.5
base = np.arange(mmin-(dm/2.),mmax+(dm/2.),dm)
base2 = base[:-1]+(base[1]-base[0])/2.

mr = U.get_data(cat,0)

m1,m2 = np.histogram(mr,base,density=False)
plt.clf()
plt.semilogy(base2,m1*total_area,'-ko')
U.put_data(outfilename,(base2,m1),'# basem mr','%.2f %.2f')
