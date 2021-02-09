__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import coeio as C
import splus_calib_tools as ct


# General root to codes
root2codes = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/'

# BPZ-related files and roots.
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
#spectra = 'COSMOSeB11new_recal'
spectra = 'eB11'

# J-PLUS related data
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/JPLUS_specz_sample_March2018/'

jplus_spz_cat = root2cats + 'jplus2sdssdr12_z05.cat'
jplus_spz_columns = root2cats + 'jplus2sdssdr12_z05.columns'
tiles = U.get_data(jplus_spz_cat,1)
single_tiles = N.unique(tiles).astype(int)
n_single_tiles = len(single_tiles)

# Reading entire catalogue.
header_jplus = C.loadheader(jplus_spz_cat)
data_jplus   = C.loaddata(jplus_spz_cat)
#C.savedata(data1[sample,:],outname,dir="",header=head1)

#n_single_tiles = 1
for ii in range(n_single_tiles):
   tile_jplus_spz_upl_cat = root2cats+'tile_%i_jplus2sdssdr12_z05.ecor.upp.cat'%(single_tiles[ii])
   good = N.equal(tiles,single_tiles[ii])
   n_gals = len(data_jplus[good,0])
   if n_gals>50:

    if not os.path.exists(tile_jplus_spz_upl_cat):
       tile_jplus_spz_cat = root2cats+'tile_%i_jplus2sdssdr12_z05.cat'%(single_tiles[ii])
       # Individual (tile) catalogue.
       if not os.path.exists(tile_jplus_spz_cat):
          good = N.equal(tiles,single_tiles[ii])
          C.savedata(data_jplus[good,:],tile_jplus_spz_cat,dir="",header=header_jplus)

       tile_jplus_spz_minerr_cat = root2cats+'tile_%i_jplus2sdssdr12_z05.ecor.cat'%(single_tiles[ii])
       if not os.path.exists(tile_jplus_spz_minerr_cat):
          ct.minimum_photouncert(tile_jplus_spz_cat,jplus_spz_columns)

       # Including upper-limits if necessary.
       ct.replace_photo_uncert(tile_jplus_spz_minerr_cat,jplus_spz_columns)

    ### ZP-recalibration
    cali_columns = tile_jplus_spz_upl_cat[:-3]+'cali.columns'
    if not os.path.exists(cali_columns):

        cmd1 = 'python %sfullcalibrator_%s.py %s '%(root2bpz,spectra,tile_jplus_spz_upl_cat)
        cmd1 += '-cols %s -outcol %s '%(jplus_spz_columns,cali_columns)
        print cmd1
        try: os.system(cmd1)
        except: continue

        # Remove temporal files created during calibration.
        cmd0 = '/bin/rm %s*temporcal*'%(root2cats)
        try: os.system(cmd0)
        except: continue

        # Save calibration plots.
        new_folder = os.path.dirname(tile_jplus_spz_upl_cat)+'/zp_plots/'
        if not os.path.exists(new_folder):
           cmd8 = '/bin/mkdir %s '%(new_folder)
           try: os.system(cmd8)
           except: continue

        nick_name_field = tile_jplus_spz_upl_cat.split('/')[-1][:-4]
        new_folder_2 = new_folder +'%s/'%(nick_name_field)
        print 'new_folder_2:',new_folder_2
        if not os.path.exists(new_folder_2):
           cmd8 = '/bin/mkdir %s '%(new_folder_2)
           try: os.system(cmd8)
           except: continue

        plots_names = root2codes+'*resoffset.png'
        cmd6 = '/bin/mv %s %s'%(plots_names,new_folder_2)
        try: os.system(cmd6)
        except: continue


    # Running BPZ
    bpz_cali  = tile_jplus_spz_upl_cat[:-3]+'cali.bpz'
    flux_cali = tile_jplus_spz_upl_cat[:-3]+'cali.flux_comparison'
    hdf5_cali = tile_jplus_spz_upl_cat[:-3]+'cali.hdf5'
    if not os.path.exists(bpz_cali):
              cmd2  = 'python %sbpz.py %s '%(root2bpz,tile_jplus_spz_upl_cat)
              cmd2 += '-COLUMNS %s -OUTPUT %s '%(cali_columns,bpz_cali)
              #cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.015 -INTERP 5 -ZMAX 0.3 '
              cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.015 -INTERP 5 -ZMAX 0.5 '
              cmd2 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(flux_cali)
              cmd2 += '-SPECTRA %s.list -ZMIN 0.0001 -PRIOR flat -DZ 0.0001 '%(spectra)
              cmd2 += '-HDF5 no '
              cmd2 += '-ZP_ERRORS "0.03,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02" '
              print cmd2
              os.system(cmd2)

