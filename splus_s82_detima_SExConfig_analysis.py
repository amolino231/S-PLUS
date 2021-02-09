__author__ = 'albertomolino'

import os,sys
import numpy as N
import useful as U
from astropy.io import fits

# Roots
SExt_config_file = 'splus.sex'

# First a couple of routines to be used.
def decapfile(filename):
    nick   = filename.split('/')[-1:][0]
    ending = nick.split('.')[-1:][0]
    return filename[:-(len(ending)+1)]

def multiply_image_bya_number(image_in,number,image_out):
    if image_in[:-2] == 'fz':
        head1 = fits.open(image_in)[0].header
        data1 = fits.open(image_in)[1].data
    else:
        head1 = fits.open(image_in)[1].header
        data1 = fits.open(image_in)[0].data
    fits.writeto(image_out,data1*number,head1,clobber=True)


# This is THE routine.
def spurious_detect_threshold(image):
    """
    It runs SExtractor twice:
    - first on the original image
    - second on the inverted image
    spanning a number of values for
    the following parameters:
     * DETECT_THRESH
     * DETECT_MINAREA
    It saves the amount of sources detected on
    each side of the image and its fraction.
    ---- USAGE:
    import splus_s82_detima_SExConfig_analysis
    from splus_s82_detima_SExConfig_analysis import *
    image = 'image.fits'
    spurious_detect_threshold(image)

    """

    # DETECT_THRESHOLD range
    det_th = N.arange(0.5,3.1,0.1)
    n_det_th = len(det_th) #Number of bins.

    # DETECT_MINAREA range
    det_mar = N.arange(1,5,1)
    n_det_mar = len(det_mar) #Number of bins.

    # Final variables where to save the results.
    fspd = N.zeros((n_det_th,n_det_mar),float)
    num_det = N.zeros((n_det_th,n_det_mar),float)

    # If necessary, it creates an inverted image.
    image_inv = decapfile(image)+'_inv.fits'
    if not os.path.exists(image_inv):
       print 'Creating an inverse image...'
       factor = -1.
       multiply_image_bya_number(image,factor,image_inv)

    #Running SExtractor
    for ii in range(n_det_th):
        for jj in range(n_det_mar):
            # New values
            threshold = det_th[ii]
            minarea = det_mar[jj]
            print '======'
            print 'Running SExtractor with:'
            print 'DETECT_THRESH: %.1f'%(threshold)
            print 'DETECT_MINAREA: %.1f'%(minarea)
            print '      '

            # It creates two new (temporal) catalogs.
            # First, using the original image
            SExt_temp_catalog_ori = decapfile(image)
            SExt_temp_catalog_ori += 'DT%.1fDA%.1f.cat'%(threshold,minarea)
            if not os.path.exists(SExt_temp_catalog_ori):
                cmd3 ='sex %s -c %s' %(image,SExt_config_file)
                cmd3 += '-CATALOG_NAME %s '%(SExt_temp_catalog_ori)
                cmd3 += '-DETECT_THRESH %.1f '%(threshold)
                cmd3 += '-DETECT_MINAREA %.1f '%(minarea)
                print cmd3
                try: os.system(cmd3)
                except: print 'Impossible to run SExtractor !!'

            # Second, using the inverted image
            SExt_temp_catalog_inv = decapfile(image_inv)
            SExt_temp_catalog_inv += 'DT%.1fDA%.1f.inv.cat'%(threshold,minarea)
            if not os.path.exists(SExt_temp_catalog_ori):
                cmd4 ='sex %s -c %s' %(image_inv,SExt_config_file)
                cmd4 += '-CATALOG_NAME %s '%(SExt_temp_catalog_inv)
                cmd4 += '-DETECT_THRESH %.1f '%(threshold)
                cmd4 += '-DETECT_MINAREA %.1f '%(minarea)
                print cmd4
                try: os.system(cmd4)
                except: print 'Impossible to run SExtractor !!'


            print 'Measuring detections...'
            if os.path.exists(SExt_temp_catalog_ori) and os.path.exists(SExt_temp_catalog_inv):
               id_o,x_o,y_o = U.get_data(SExt_temp_catalog_ori,(0,1,2))
               id_i = U.get_data(SExt_temp_catalog_inv,0)
               # Excluding image edges to avoid other effects.
               good  = N.greater(x_o,1500.) * N.less(x_o,10300.)
               good *= N.greater(y_o,1100.) * N.less(y_o,10000.)
               id_o_redu = N.compress(good,id_o)
               id_i_redu = N.compress(good,id_i)
               # Ratio among detections on each side.
               fspd[ii,jj] = 100.*(len(id_i_redu)/(1.*len(id_o_redu)))
               num_det[ii,jj] = len(id_o_redu)
            else:
                print 'Catalogs do not exist!!'
                sys.exit()

    ## Saving final results into an external file.
    # Fraction of spurious detections
    out_filename_1 = decapfile(image)+'f.mat'
    U.put_2Darray(out_filename_1,fspd)
    # Number of detections on original image.
    out_filename_2 = decapfile(image)+'n.mat'
    U.put_2Darray(out_filename_2,num_det)


