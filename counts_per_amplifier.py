__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import matplotlib.pyplot as plt

"""
This routines looks for statistical deviation in the number counts
of detected sources in the different amplifiers. It uses a list of
catalogues as input and returns a plot and a ascii file.
"""


mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
final_root = mainroot+'data_quality/'
root2cats = mainroot + 'splus_cats/'
lista_cats = 'photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)

out_filename = final_root+'countsxamplif'

if not os.path.exists(out_filename):
   # Variables where to store the information
   counts = N.zeros((2,8),float) # 2x8 amplifiers.

   for ii in range(n_cats):
       catalog = cats_names[ii]
       print 'Reading file: ',catalog.split('/')[-1]
       x,y,mr = U.get_data(catalog,(3,4,78))
       good_sample  = N.greater_equal(mr,14)
       good_sample *= N.less_equal(mr,20)
       x,y,mr = U.multicompress(good_sample,(x,y,mr))

       #Defining edges
       min_x = min(x)
       max_x = max(x)
       min_y = min(y)
       max_y = max(y)
       dx = (max_x-min_x)/8.
       dy = (max_y-min_y)/2.

       #print 'dx: ',dx
       #print 'dy: ',dy

       # Starting loop
       for jj in range(2):
           for hh in range(8):
               if hh<1:
                   x_sel = N.less_equal(x,min_x+(hh+1)*dx)
               else:
                   x_sel  = N.greater_equal(x,min_x+(hh*dx))
                   x_sel *= N.less_equal(x,min_x+(hh+1)*dx)
               if jj<1:
                  y_sel = N.less_equal(y,min_y+dy)
               else:
                  y_sel = N.greater_equal(y,min_y+dy)

               # Counting sources.
               counts[jj,hh] += len(x[x_sel*y_sel])

   max_value = counts.max()

   cc = N.reshape(counts/(1.*max_value),(1,16))
   base = N.arange(16)+1
   #Saving data.
   U.put_data(out_filename+'.txt',(base,cc[0]/(1.*n_cats)),'# amp norm_counts')
else:
    base,cc = U.get_data(out_filename+'.txt',(0,1))


# Starting plot.
plt.figure(2, figsize=(10,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
plt.bar(base-0.2,cc[0]/(1.),0.4,color='grey',alpha=0.8,linewidth=2)
plt.xlabel('Amplifier',size=24)
plt.ylabel('Relative Number Counts',size=26,labelpad=6)
plt.xticks(fontsize=22)
plt.yticks(fontsize=20)
plt.xlim(0.1,16.9)
plt.ylim(0.9,1.0)
plt.grid()

plt.savefig(out_filename+'.png',dpi=90)




