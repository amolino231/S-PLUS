__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U

"""
This routine is executed to run BPZ on its OT=yes mode
and obtained the expected vs observed fluxes for all galaxies
with spec-z info. These files will be used then to compute re-calibrations
to the initial SED library using 'splus_SEDrecal.py'.
Finally, new ZP-corrections will be performed for the new SED library
and BPZ will run again over all catalogues.

"""


mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/'
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root2codes = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/'
lista_cats = 'photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)
apertures = ['auto']
n_apers = len(apertures)

### spectral library
spectra = 'COSMOSeB11new_recal'

final_root = root2cats+'%s_OT_after_SEDrecal/'%(spectra) # New root
if not os.path.exists(final_root):
    cmd8 = '/bin/mkdir %s '%(final_root)
    os.system(cmd8)

root_to_columns = root2cats+'%s/'%(spectra)

# splus_columns = root2cats + 'splus11_auto.columns'
sdss_s82_spz_cat = root2cats + 'S82_SDSS_z05.cat'


########## Cross-matching with SDSS/S82.
for ggg in range(n_cats):
    splus_speczcat = cats_names[ggg][:-3] + 'spz.z05.cat'
    for hhh in range(n_apers):
        nickname = os.path.basename(splus_speczcat)
        cali_cols = root_to_columns + nickname[:-3]+'%s_cali.columns'%(apertures[hhh])
        bpz_cali  = final_root+nickname[:-3]+'%s_cali.OT.bpz'%(apertures[hhh])
        flux_cali = final_root+nickname[:-3]+'%s_cali.OT.flux_comparison'%(apertures[hhh])
        #hdf5_cali = final_root+nickname[:-3]+'%s_cali.OT.hdf5'%(apertures[hhh])
        if not os.path.exists(bpz_cali):
            cmd2  = 'python %sbpz.py %s '%(root2bpz,splus_speczcat)
            cmd2 += '-COLUMNS %s -OUTPUT %s '%(cali_cols,bpz_cali)
            cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.015 -INTERP 5 -ZMAX 0.5 '
            cmd2 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(flux_cali)
            cmd2 += '-SPECTRA %s.list -ZMIN 0.001 -PRIOR flat -DZ 0.0001 '%(spectra)
            cmd2 += '-HDF5 no -ONLY_TYPE yes '
            cmd2 += '-ZP_ERRORS "0.03,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02" '
            os.system(cmd2)


        """
        # Renaming HDF5 files
        if not os.path.exists(hdf5_cali):
           temporal_hdf5 = os.path.dirname(splus_speczcat)+'/'
           temporal_hdf5 += os.path.basename(splus_speczcat).split('.')[0] + '.hdf5'
           print 'temporal_hdf5',temporal_hdf5
           print 'hdf5_cali',hdf5_cali
           cmd7 = '/bin/mv %s %s '%(temporal_hdf5,hdf5_cali)
           print cmd7
           os.system(cmd7)
        """

        # Remove temporal files created by BPZ.
        cmd3 = '/bin/rm %s*.npy'%(root2cats)
        try: os.system(cmd3)
        except: continue



