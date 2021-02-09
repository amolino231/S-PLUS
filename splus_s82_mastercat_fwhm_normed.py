__author__ = 'albertomolino'

import os,sys
import numpy as N
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U
import coeio as C
import splus_calib_tools as sct


# Roots & catalog.
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cat = mainroot + 'released_catalogues/'
#mastercat = root2cat + 'SPLUS_STRIPE82_master_catalogue_edr_march2018.cat'
mastercat = root2cat + 'S82_master_gals.cat'

"""
# This loop was needed for the original version
# where the FIELDS were strings rather than integers.
tiles_name = U.get_str(mastercat,0)
nt = len(tiles_name)
tiles = N.zeros(nt)
for ss in range(nt):
    tiles[ss] = int(tiles_name[ss].split('-')[-1])
"""
tiles = U.get_data(mastercat,0) # use get_str with original cat!
single_tiles = N.unique(tiles).astype(int)
n_single_tiles = len(single_tiles)

# Reading entire catalogue.
header_mastercat = C.loadheader(mastercat)
data_mastercat   = C.loaddata(mastercat)

for ii in range(n_single_tiles):
    tile_cat = root2cat+'tile%i_S82.cat'%(single_tiles[ii])
    good = N.equal(tiles,single_tiles[ii])
    n_gals = len(data_mastercat[good,1])
    if not os.path.exists(tile_cat):
       """
       fwhm = data_mastercat[good,9]
       magr = data_mastercat[good,82] # Petro
       seeing,stars = sct.get_seeing_from_data_pro(fwhm,magr)
       fwhm_norm = fwhm /(1.*seeing)
       data_mastercat[good,9] = fwhm_norm
       data_mastercat[good,0] = N.ones(n_gals) * int(single_tiles[ii])
       """
       C.savedata(data_mastercat[good,:],tile_cat,dir="",header=header_mastercat)
