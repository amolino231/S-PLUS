__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import jplus_calib_tools as jct
#import alhambra_photools
#from alhambra_photools import reset_zperrors

"""
1. Open the JPLUS/SDSS catalogue
2. Separate the catalogue in Tiles
3. If #>ng_min: calibrate
4. Run BPZCAL not changing FixFilters
5. Run BPZ w/ & w/o ZP calibration.
6. Save BPZfiles
"""

# Paths definition
root2maincat = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Sept2017/'
maincat = root2maincat+'SV02_March07.lite.spz.cat'
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
tile_root = root2maincat+'tiles/'
zpe = 0.03 # minimum zp_error for cal.columns
if not os.path.exists(tile_root):
   os.system('/bin/mkdir '+tile_root)


nocal_tiles = tile_root+'bpzcal.info.txt'
tiles_info = open(nocal_tiles,'w')
tiles_info.write('# Tiles Ng seeing \n')

# catalogues parameters
header_cat = '# 1.tile_id 2.object_id 3.ra 4.dec 5.fwhm 6.uJAVA 7.uJAVA_err 8.gSDSS \
 9.gSDSS_err 10.rSDSS 11.rSDSS_err 12.iSDSS 13.iSDSS_err 14.zSDSS 15.zSDSS_err \
 16.J0395 17.J0395_err 18.J0410 19.J0410_err 20.J0430 21.J0430_err \
 17.J0515 18.J0515_err 19.J0660 20.J0660_err 21.J0861 22.J0861_err \
 23.zs_SDSS '

# BPZCAL parameters
columns11 = root2maincat+'jplus11.columns' #ID==Tile
columns11b = root2maincat+'jplus11b.columns' #ID==ID
ngal_min = 25
#FixFilters = ["jplus_r_sdss.res"]

# Opening main.cat and generate tile.cat
tile_num,obj_id,ra,dec,fw = U.get_data(maincat,(0,1,2,3,4))
u_m, u_em, g_m, g_em, r_m, r_em = U.get_data(maincat,(5,6,7,8,9,10))
i_m, i_em, z_m, z_em, f395_m, f395_em = U.get_data(maincat,(11,12,13,14,15,16))
f410_m, f410_em, f430_m, f430_em, f515_m, f515_em = U.get_data(maincat,(17,18,19,20,21,22))
f660_m, f660_em, f861_m, f861_em, zpec = U.get_data(maincat,(23,24,25,26,28))

only_tiles = jct.num_tiles(tile_num) # tiles w/o duplications
n_tiles = len(only_tiles) # number of Tiles
print 'Total number of Tiles to calibrate: ',n_tiles

for ii in range(n_tiles):
    ref_tile = only_tiles[ii]
    gs = N.less(abs(tile_num-ref_tile),1)# GoodSample
    tile_cat = tile_root+'%i.cat'%(ref_tile)
    U.put_data(tile_cat,
               (tile_num[gs],obj_id[gs],ra[gs],dec[gs],fw[gs],
            u_m[gs], u_em[gs], g_m[gs], g_em[gs], r_m[gs], r_em[gs],
            i_m[gs], i_em[gs], z_m[gs], z_em[gs], f395_m[gs], f395_em[gs],
            f410_m[gs], f410_em[gs], f430_m[gs], f430_em[gs],
            f515_m[gs], f515_em[gs], f660_m[gs], f660_em[gs],
            f861_m[gs], f861_em[gs], zpec[gs]),
               header_cat)

    if os.path.exists(tile_cat):
       ids = U.get_data(tile_cat,0)
       ngals = len(ids)
       seeing = jct.get_seeing_from_data(fw[gs],r_m[gs]) #Estimates the seeing
       tiles_info.write('%s %i %.3f \n'%(tile_cat,ngals,seeing))
       if ngals>=ngal_min:
          cali_columns = tile_cat[:-3]+'cal.columns'
          if not os.path.exists(cali_columns):
             cmd1 = 'python %sfullcalibrator_amb.py %s '%(root2bpz,tile_cat)
             cmd1 += '-cols %s -outcol %s '%(columns11b,cali_columns)
             print cmd1
             os.system(cmd1)
          else:
              print 'The file already exists!'
       else:
           print 'Few galaxies to calibrate'

    # here we run BPZ on Tile.cat using default COLUMNS.
    bpz_tile_cat = tile_cat[:-3]+'nocal.bpz'
    if not os.path.exists(bpz_tile_cat):
       bpz_flux_tile_cat = tile_cat[:-3]+'flux_comparison'
       cmd2  = 'python %sbpz.py %s '%(root2bpz,tile_cat)
       cmd2 += '-COLUMNS %s -OUTPUT %s '%(columns11b,bpz_tile_cat)
       cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.03 -INTERP 5 -ZMAX 1.0 '
       cmd2 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(bpz_flux_tile_cat)
       cmd2 += '-SPECTRA GOSMOSeB11.list -PRIOR flat -ZMIN 0.001'
       os.system(cmd2)

    # here we run BPZ on Tile.cat using calibrated COLUMNS.
    if os.path.exists(cali_columns):
        # bpz_tile_cal_cat = tile_cat[:-3]+'cal.bpz'
        bpz_tile_cal_cat = tile_cat[:-3]+'GOSMOSeB11.cal.bpz'
        if os.path.exists(cali_columns):
           bpz_flux_tile_cat = tile_cat[:-3]+'cal.flux_comparison'
           cmd4  = 'python %sbpz.py %s '%(root2bpz,tile_cat)
           cmd4 += '-COLUMNS %s -OUTPUT %s '%(cali_columns,bpz_tile_cal_cat)
           cmd4 += '-CHECK yes -SIGMA_EXPECTED 0.03 -INTERP 5 -ZMAX 1.0 '
           cmd4 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(bpz_flux_tile_cat)
           cmd4 += '-SPECTRA GOSMOSeB11.list -PRIOR flat -ZMIN 0.001 '
           cmd4 += '-ZP_ERRORS "%.3f,%.3f,%.3f,%.3f,%.3f,'%(zpe,zpe,zpe,zpe,zpe)
           cmd4 += '%.3f,%.3f,%.3f,%.3f,%.3f,'%(zpe,zpe,zpe,zpe,zpe)
           #cmd4 += '%.3f,%.3f"'%(zpe,zpe)
           cmd4 += '%.3f"'%(zpe) #Only 11 filters this run.
           print cmd4
           os.system(cmd4)


    # Removing temporal files...
    cmd3 = '/bin/rm -rf %s/*temporcal*'%(tile_root)
    #os.system(cmd3)
    cmd31 = '/bin/rm -rf *resoffset.png'
    os.system(cmd31)
    cmd32 = '/bin/rm -rf %s/*npy'%(tile_root)
    os.system(cmd32)
    print 'Deleting temporal files...'


tiles_info.close()
