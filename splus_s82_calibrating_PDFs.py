__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U

# roots
root2cats  = '/Volumes/CLASH/S82/specz/'
root2bpz   = '/Users/albertomolino/codigos/bpz-1.99.2/'
spectra    = 'COSMOSeB11new_recal'

# Reading lists
bpz_cats_list = 'bpz_short.list'
bpz_cats = U.get_str(root2cats + bpz_cats_list,0)
n_cats = len(bpz_cats)

# SIGMA_VALUES to be tested
sigma_values = ['015','020','025','030','001','10']
n_sigma_values = len(sigma_values)

for ii in range(n_sigma_values):
    final_root = root2cats + 'SC/'
    final_root += 'sc%s/'%(sigma_values[ii])
    print final_root
    if not os.path.exists(final_root):
       cmd8 = '/bin/mkdir %s '%(final_root)
       os.system(cmd8)

    for ss in range(n_cats):
        hdf5_file = final_root+os.path.basename(bpz_cats[ss])[:-4]
        hdf5_file += '.sc%s.hdf5'%(sigma_values[ii])
        print 'hdf5_file',hdf5_file
        if not os.path.exists(hdf5_file):
           bpz_file = hdf5_file[:-4]+'bpz'

           if ii<1:
               columns  = root2cats+os.path.basename(hdf5_file[:-10])+'columns'
               cat_file = root2cats+os.path.basename(hdf5_file[:-15])+'cat'
           else:
               columns  = root2cats+os.path.basename(hdf5_file[:-9])+'columns'
               cat_file = root2cats+os.path.basename(hdf5_file[:-14])+'cat'

           cmd2  = 'python %sbpzConvolved.py %s '%(root2bpz,cat_file)
           cmd2 += '-COLUMNS %s -OUTPUT %s -CONVOLVE_P yes '%(columns,bpz_file)
           cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.015 -INTERP 5 -ZMAX 1.0 '
           #cmd2 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(flux_cali)
           cmd2 += '-SPECTRA %s.list -ZMIN 0.005 -PRIOR SM -DZ 0.001 '%(spectra)
           cmd2 += '-ZP_ERRORS "0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02" '
           cmd2 += '-HDF5 yes -SIGMA_CONVOLVE 0.%s '%(sigma_values[ii])
           print cmd2
           os.system(cmd2)
           #pausa = raw_input('paused')

           temporal_hdf5 = final_root+os.path.basename(hdf5_file[:-14])+'hdf5'
           print 'temporal_hdf5',temporal_hdf5
           if os.path.exists(temporal_hdf5):
              cmd3 = '/bin/mv %s %s '%(temporal_hdf5,hdf5_file)
              os.system(cmd3)


