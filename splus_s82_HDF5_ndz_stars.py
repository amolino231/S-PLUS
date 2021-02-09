__author__ = 'albertomolino'

"""
This routine is essentially identical to 'splus_s82_HDF5_ndz.py' but
it only selects potential stars! (not galaxies).
The idea is to understand where in z stars tend to clusters, identifying
potential redshift intervals where the contamination could be higher/lower.
--
Important for cosmology or statistical analysis.

"""

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import splus_s82_hdf5_tools as to

# General roots.
root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
root += 'S82/Dec2017/'

# Reading HDF5 files.
root_to_hdf5 = root + 'splus_cats_NGSL/*/' # TBC!!
final_root_hdf5 = root_to_hdf5 + 'PDFs/'
if not os.path.exists(final_root_hdf5):
   cmd = '/bin/mkdir %s '%(final_root_hdf5)
   os.system(cmd)
hdf5list = root_to_hdf5 + 'hdf5.list'
hdf5_files = U.get_str(hdf5list,0)
n_hdf5 = len(hdf5_files)

# Reading photometric catalogues.
root_to_cats = root + 'released_catalogues/'
cat_list = root_to_cats + 'master_SPLUS_STRIPE82_photo_BPZ.list'
global_cats = U.get_str(cat_list,0)
n_cats = len(global_cats)

# Checking dimensionality!
if n_cats != n_hdf5:
   print 'Dimension missmatch!'
   sys.exit()

# Magnitude bins.
mag_bins = [14,16,18,20]
n_mags = len(mag_bins)

# Starting loop.
for ii in range(n_mags):
    print 'Processing detections with R<%i: '%(mag_bins[ii])
    final_pdf_file = final_root_hdf5 + 'master_R%i_PDF.stars.txt'%(mag_bins[ii])
    if not os.path.exists(final_pdf_file):
       for ss in range(n_cats):
           print 'reading file %i/%i '%(ss+1,n_cats)
           #1. Select potential stars with R<Ri
           catalog = global_cats[ss]
           p_star = U.get_data(catalog,133) # Probability of being a star.
           # Reading P(z) from HDF5 file.
           zz,pr,pb,pg = to.getPDF_by_mag_and_weights(hdf5_files[ss],mag_bins[ii],p_star)
           # Storing P(z)
           if ss<1:
              base_z = zz
              final_pdf_red = pr
              final_pdf_blue = pb
              final_pdf_global = pg
           else:
              final_pdf_red += pr
              final_pdf_blue += pb
              final_pdf_global += pg

       # Saving P(z|R<Ri)
       U.put_data(final_pdf_file,(zz,final_pdf_red,final_pdf_blue,final_pdf_global),'# z Pr Pb Pg ')


### Here we can do some plots.



