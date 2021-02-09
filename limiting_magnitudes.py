#! /usr/local/bin python
# -*- coding: iso-8859-1 -*-

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
from astropy.io import fits
import matplotlib.pyplot as plt

filtros = ['uJAVA','J0378','J0395','J0410','J0430','gSDSS',
           'J0515','rSDSS','J0660','iSDSS','J0861','zSDSS']


def check_rms_SPLUS(segmentation,photometry,minrad,maxrad,totnum,plots,verbose):

    """
    It describes the photometric noise in images by launching apertures
    in blank areas and describing the area vs rms dependency.

    Philosophy:
    Several apertures (of random radius and positions) will be created (over blank areas) on 
    a 'segmentation' image to estimate the real photometric error of an input 'photometry' image. 
    --------
    segmentation: Segmentation-like image (SExtractor output) used to select 'blank' areas.
    photometry:   scientific image over which estimate the real photometric errors      
    minrad,maxrad = minimun & maximum radius of the used apertures (pixels).
    totnum = total number of apertures.
    area,rms = final outputs.
    ---------
    USAGE:

---------
import script_dcluster_tools as to
segmentation = 'f814.seg.fits'
photometry   = 'f814.fits'
apertures,finalbackg,finalmeans,fluxes = to.check_rms_JPLUS(segmentation,photometry,1,21,5.0e+04,'yes',False)
----

    """ 
    if not os.path.exists(segmentation):
       print
       print 'Image %s does not exist!' %(segmentation)
       sys.exist()

    if not os.path.exists(photometry):
       print
       print 'Image %s does not exist!' %(photometry)
       sys.exist()
    
    if verbose==True: verba=1
    else: verba=0	       

    # Reading data from images
    if photometry[:-2]=='fz':
       photima = fits.open(photometry)[1].data
    else:
       photima = fits.open(photometry)[0].data
    
    if segmentation[:-2]=='fz':
       segima  = fits.open(segmentation)[1].data  
    else:
       segima  = fits.open(segmentation)[0].data

    # Final root where to save the data
    final_path = os.path.os.path.dirname(photometry)
    base_name = os.path.os.path.basename(photometry)
    len_extension = len(base_name.split('.')[-1])+1
    file_root = final_path+'/Apertures/%s'%(base_name[:-len_extension])
       
    # Physical limits (pixel) for the segmentation image.
    xGC = N.shape(segima)[1]/2.
    yGC = N.shape(segima)[0]/2.
    # For CLASH, due to the rotating frames, a maximum radius is set.
    min_maxim_radius = min([xGC,yGC])
    radialmax = min_maxim_radius # Maximum radial distance to the center [pixels]
    
    # Here the random position (X,Y) are limited in range.
    minpix_X = 1500 
    maxpix_X = N.shape(segima)[1]-1500 
    minpix_Y = 1500 
    maxpix_Y = N.shape(segima)[0]-1500 
    binhisto = 100  
    # Final vector with positions.
    x_values = N.arange(minpix_X,maxpix_X,1)
    y_values = N.arange(minpix_Y,maxpix_Y,1)
    # Total dimension for the input variables.
    maxdim = int(20*(10.+(4*(((maxrad)*(maxrad+1))/2.))))
    # Defining other variables.
    XX = N.zeros((maxdim),float)
    YY = N.zeros((maxdim),float)
    XO = N.zeros((maxdim),float)
    YO = N.zeros((maxdim),float)
    RR = N.zeros(totnum)
    # Range of apertures to be launched.
    apertures = N.arange(minrad,maxrad,1)
    n_apertures = len(apertures)
    # Length definition for the final outputs.
    finalbackg = N.zeros(n_apertures,'float64')
    finalmeans = N.zeros(n_apertures,'float64') 
    gausshisto = N.zeros((binhisto,2*n_apertures),dtype='float64') 
    # Starting the analysis.
    mmm = 0
    for app in range(n_apertures):
        raper = apertures[app] 
        # minrad = maxrad = raper
        if verba:
           print 'Interation %i out of %i ' %(app+1,n_apertures)
           print '-------------------------'
            
        # New temporal variables (erased in every loop).
        fluxes = N.zeros(totnum,dtype='float64')
        sbackg = N.zeros(totnum,dtype='float64') 
        ff = -1
        
        # Now it runs until it gets "totnum" measurements.
        hh = -1
        contador = 0
        while contador < (totnum):
            kk = -1
            # Random x,y numbers to estimate the position
            # to place the aperture.
            xo = N.random.random_integers(minpix_X,maxpix_X)  
            yo = N.random.random_integers(minpix_Y,maxpix_Y)
            if verba: print 'xo,yo',xo,yo
            # Corresponding radial distance to the center.
            tempradii = N.sqrt((xo-xGC)*(xo-xGC)+(yo-yGC)*(yo-yGC))
            if tempradii < (radialmax+1):
               # Now it is computed the shape of the aperture. 
               Xr = N.zeros((raper*raper),float)
               Yr = N.zeros((raper*raper),float)
               for ii in range(raper):
                   xvalue = xo+ii
                   for jj in range(raper):
                       kk += 1
                       yvalue = yo+jj
                       Xr[kk] = xvalue
                       Yr[kk] = yvalue
                       
               # Here it checks the blanckness of the aperture.                                 
               tempflux = area2noise(segima,photima,Xr,Yr)
               if raper<2 :
                  if tempflux != -999. : 
                     fluxes[hh] = tempflux
                     if verba: print 'Adding flux: ',tempflux
                     hh += 1 
                     contador += 1
                         
               if raper>1:
                  if tempflux[0] != -999. :
                     fluxes[hh] = tempflux.sum() - (finalmeans[0] * raper**2) #why?
                     contador += 1
                     hh += 1 
                     if verba: print 'hh',hh
                   
        # Computing values from the sample.    
        sigfluxes = U.std_robust(fluxes) 
        good      = U.less_equal(abs(fluxes),5.*sigfluxes)
        fluxes    = U.compress(good,fluxes)
        
        # Storing the background dispersion & mean inside that aperture.   
        finalbackg[app] = U.std(fluxes)
        finalmeans[app] = U.mean(fluxes)
        
        if plots == 'yes': 
           plt.figure(1,figsize = (7,6),dpi=70, facecolor='w', edgecolor='k')
           plt.clf()
           # va1,va2 = N.histogram(fluxes,binhisto,normed=1)
           va1,va2,va3 = plt.hist(fluxes,binhisto,normed=1,facecolor='black',alpha=0.5,linewidth=1.5)
           baseh = va2[0:-1] + ((va2[1]-va2[0])/2.)
           nele  = len(fluxes)
           mu    = U.mean(fluxes)
           sig   = U.std(fluxes)
           # yh    = U.normpdf(va2,mu,sig)
           # plt.plot(va2,yh,'r-',linewidth=3,alpha=0.7)
           mu = U.mean_robust(fluxes) # repeated
           sig = U.std_robust(fluxes) # repeated
           # yh = U.normpdf(va2,mu,sig) # repeated
           # plt.plot(va2,yh,'r--',linewidth=3,alpha=0.7) # repeated          
           plt.legend([('MEAN: %.4f ''\n'' RMS:  %.4f '%(mu,sig)),
                       'Aperture: %i $pix$'%(raper*raper)],
                       numpoints=1,loc='upper right',fontsize=14)
           plt.xlim(mu-4*sig,mu+4*sig)
           plt.xlabel('Aperture Flux [ADU]',size=20)
           plt.ylabel('Number Counts',size=20)
           plt.xticks(fontsize=17),plt.yticks(fontsize=17)
           nameima = photometry.split('/')[-1:][0]
           plt.ylim()
           figure2name = file_root+'_hfaper_%i.png' %(raper)
           plt.savefig(figure2name,dpi=150)
           plt.close()
           
           # Here it saves the info from the histogram.
           ind1 = app*2
           ind2 = app*2+1
           if verba: print 'ind1,ind2',ind1,ind2
           gausshisto[:,ind1] = baseh # va2
           gausshisto[:,ind2] = va1   # yh
            
    # At this point all apertures have been computed.        
    # Now it will represent the sigma_vs_area dependency.
    sigmas = finalbackg #-abs(finalmeans)
    aa,bb = sigmafit(sigmas,sigmas[0],apertures)
      
    if plots == 'yes':
       plt.figure(2, figsize = (7,6),dpi=80, facecolor='w', edgecolor='k')
       plt.clf()
       plt.plot(apertures,sigmas[0]*apertures*(aa+bb*apertures),'k-',apertures,apertures*sigmas[0],'r-')
       plt.legend([('%.3f$\sqrt{N}$ (%.3f + %.3f$\sqrt{N}$)' %(sigmas[0],aa,bb)),
                   '%.3f$\sqrt{N}$ | Poisson Distribution '%(sigmas[0])],
                   numpoints=1,loc='upper left')
       plt.plot(apertures,sigmas,'ko')
       plt.xlim(0.,max(apertures)+1)
       plt.xlabel('$\sqrt{N}$',size=18)
       plt.ylabel('$\sigma$',size=20)
       plt.xticks(fontsize=15)
       plt.yticks(fontsize=15)
       nick2 = photometry.split('/')[-1:][0]

       figure1name = file_root+'_apersigma.png'
       plt.savefig(figure1name,dpi=150)
       plt.close()

    # Saving outputs in ASCII files.
    fileout = file_root+'.apertures.txt'
    header = '# AREA[pix] RMS(std[counts]) MEAN(mean_robust[counts])'
    U.put_data(fileout,(apertures,sigmas,finalmeans),header)
    #print 'Saving data... in %s' %(fileout)

    """
    gausstxt = photometry[:-5]+'_gaussfit.txt'
    outfile = open(gausstxt,'w')
    kk = 0
    header = ['gausshisto_X','gausshisto_Y']  
    for ss in range(n_apertures):
        for hh in range(len(header)):
            kk +=1
            outfile.write('%s' ' %i ' ' %s_%i ' ' \n' %('# ', kk,header[hh],ss+1))
        outfile.write('%s' ' \n' %('#'))
        
    for ii in range(binhisto):
        for jj in range(2*n_apertures):  
            outfile.write('%.4f ''\t' %(gausshisto[ii,jj]))
        outfile.write('\n') 
    outfile.close()
    """
    
    #return apertures,finalbackg,finalmeans,fluxes




def area2noise(segima,photima,Xpix,Ypix):
    """
    It measures the fluxes, inside the aperture, in the case
    those pixels have no signal on the sementation image.
    --------------------------------------------------------
    segima & photima are not strings but the matrix itself !!

    """
    nx = len(Xpix)   # Number of X-pixels inside each circle.
    ny = len(Ypix)   # Number of Y-pixels inside each circle.
    flux = N.zeros(nx) 
    for efe in range(1): 
        raw = 0
        if nx == 1 :
           # Matrix elements
           hh = int(Xpix[0])-1
           gg = int(Ypix[0])-1
           # Assoc. values inside matrices
           # print 'segima[gg,hh]',segima[gg,hh]
           raw += int(segima[gg,hh])
           flux = float(photima[gg,hh])

           if raw == 0 :  
              #print 'Good region. Measuring background...'  
              fluxes = flux
           else:
              fluxes = - 999. 
    
        else:
           for hhh in range(nx):
                hh = int(Xpix[hhh])-1
                gg = int(Ypix[hhh])-1                
                # Assoc. values inside matrices
                raw += float(segima[gg,hh])
                flux[hhh] = float(photima[gg,hh])
           if raw == 0 :  
              #print 'Good region. Measuring background...'  
              fluxes = flux
           else:
              fluxes = (flux * 0.) - 999. 
                
    return fluxes     



def sigmafit(sigma,sigma1,sqarea):

    """
    It linearily fits the input data as follows:
    y = bx + a, where:   y = sigma/(sigma1*sqarea) 
                         x = sqarea
    ---------------------------------------------------- 
    All three inputs sigma, sigma1 & sqarea are 1D vectors.
    """
    # Definition of the linear equation
    print N.shape(sigma),N.shape(sigma1),N.shape(sqarea)
    sigma1v = N.ones(len(sigma))*sigma1
    y = sigma / (sigma1v * sqarea)
    x = sqarea
    # Fitting data linearily.
    ccx = U.lsq(x,y)
    bb = ccx.b      # Coefficient prop. to N*N
    aa = ccx.a      # Coefficient prop. to N 

    return aa,bb





def combineimages(lista,outfile,normed):

    """
    It combines a list of image using Fits.
    The final header information will be the same as
    the image selected (1 by default).
    
    """
    # Declaring some variables.
    print 'Number of images to be combined: %i' %(len(lista))
    print 'Starting it out...' 

    headim = 1
    head = fits.open(lista[0])[0].header
    base = fits.open(lista[0])[0].data
    dimx = len(base[0,:])
    dimy = len(base[:,0])
    matrix = N.zeros((dimx,dimy),float)  # New (final) matrix. 

    for ss in range(len(lista)):
        temp  = fits.open(lista[ss])[0].data
        tempx = len(temp[0,:])
        tempy = len(temp[:,0])

        if (tempx == dimx) and (tempy == dimy):
           matrix += temp
        else:
            print 'Image %i has a different dimension!' %(ss+1)
            sys.exit()   

    # If normed, it divides the final matrix
    # by the number of combined images.
    # Useful for weight-maps.
    if normed:
       matrix /= float(len(lista)) 
            
    # Now it saves the final image
    fits.writeto(outfile,matrix,head)


def JPLUS_image_gain(image):
    """
    It reads the image's header and
    look for the parameter
    which accounts for the gain value.
    """
    head = fits.open(image)[0].header   
    return float(head['GAIN'])       


def get_JPLUS_gains(scimas):
    """
    it gets the gain values for a list of JPLUS images.
    """
    imas  = U.get_str(scimas,0)
    nimas = len(imas)
    vals = N.zeros(nimas)
    for ss in range(nimas):
        vals[ss]=JPLUS_image_gain(imas[ss])
    return vals    

def decapfile(filename):
    """
    It removes the file extension.
-----
file = '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg.fits'
decapfile(file) ==> '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg'
    """
    nick   = filename.split('/')[-1:][0] 
    ending = nick.split('.')[-1:][0]
    dim    = len(ending)+1
    name   = filename[:-dim]
    return name


def get_nickname(file):

    """
    It purges both the root and the extension from a file
-----
file = '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg.fits'
get_nickname(file) ==> 'f04p01_1_acs.deg'
    """
    pepe = file.split('/')[-1:][0].split('.')[:-1][0]
    return pepe



def get_list_detect_imas(cluster):
    """
    Here i need to include more info to
    separate weight than drz images.
    """
    root2images =finalroot+'images/%s/'%(cluster)
    # R-band
    cmd1  = '/bin/ls %s*rSDSS*swp.fits '%(root2images)
    cmd1 += '> %sdetima.riz.temp'%(root2images)
    os.system(cmd1)
    # I-band
    cmd2  = '/bin/ls %s*iSDSS*swp.fits '%(root2images)
    cmd2 += '>> %sdetima.riz.temp'%(root2images)
    os.system(cmd2)    
    # z-band
    cmd3  = '/bin/ls %s*zSDSS*swp.fits '%(root2images)
    cmd3 += '>> %sdetima.riz.temp'%(root2images)
    os.system(cmd3)          
    # Reading final list
    finalist = '%sdetima.riz.temp'%(root2images) 
    if not os.path.exists(finalist): 
       print 'List of detection images not created!'
       sys.exit()
    else:
        dimas = U.get_str(finalist,0)
    return dimas


def get_list_detect_wimas(cluster):
    """
    Here i need to include more info to
    separate weight than drz images.
    """
    root2images =finalroot+'images/%s/'%(cluster)
    # R-band
    cmd1  = '/bin/ls %s*rSDSS*weight*fits '%(root2images)
    cmd1 += '> %sdetwima.riz.temp'%(root2images)
    os.system(cmd1)
    # I-band
    cmd2  = '/bin/ls %s*iSDSS*weight*fits '%(root2images)
    cmd2 += '>> %sdetwima.riz.temp'%(root2images)
    os.system(cmd2)    
    # z-band
    cmd3  = '/bin/ls %s*zSDSS*weight*fits '%(root2images)
    cmd3 += '>> %sdetwima.riz.temp'%(root2images)
    os.system(cmd3)          
    # Reading final list
    finalist = '%sdetwima.riz.temp'%(root2images) 
    if not os.path.exists(finalist): 
       print 'List of detection images not created!'
       sys.exit()
    else:
        wimas = U.get_str(finalist,0)
    return wimas


def correct_SExt_uncertainties(cluster):
# (catalog,columns,zpts,gains,area2rms,weightimas,arinarout,finalcat):
    """
    It reads the input catalogue and corrects*
    its photometric errors empirically using direct area_v_sigma estimations*.
    It returns the name of the new and corrected catalog. 
    --
    *** THE POSITION OF AREA,MAGS & ERRMAGS NEED TO BE CHECK OUT BEFORE RUNNING IT.
    ----------------------------------------------------------------------------
    REQUIREMENTS:
      - An input catalogue (catalog)
      - Its corresponding COLUMNS file (columns)
      - A list with all ZPS (2nd column) (zpts)
      - A list with all GAINS (2nd column) (gains)
      - A list with all area2sigma files (1 per band) (area2rms)
      - A list with all WEIGHT-MAPS (weightimas)
    """
    # If weight = 1, it uses the Weight-maps to calculate photo.uncertainties.
    weight = 0
    # If verbose = 1, additional information is displayed during the analysis.
    verbose = 0
    # This factor serves to plot several check figures.
    check1 = 1
    check2 = 1

    catalog = finalroot +'/%s/images/%s.photo.cat'%(cluster,cluster)
    columns = finalroot +'/%s/images/%s.photo.columns'%(cluster,cluster)
    scimas  = finalroot+'/%s/images/sci.list'%(cluster)
    weightimas = finalroot+'/%s/images/wht.list'%(cluster)
    zpts = U.get_data(finalroot+'/%s/images/%s.zpt.cat'%(cluster,cluster),0)
    gains = get_JPLUS_gains(scimas)

    # Name for the final corrected catalogue
    newphotcat = finalroot +'/%s/images/%s.photo.err.cat'%(cluster,cluster)

    # To account for differences in exposure time, we need a list with ALL WEIGHT-maps.
    if weight: wimas = U.get_str(weightimas,0)

    if verbose: print 'Catalog: ',catalog
    if os.path.exists(catalog):
       mm  = get_magnitudes(catalog,columns)
       em = get_errmagnitudes(catalog,columns)
       xx,yy,aper = U.get_data(catalog,(3,4,5))
       data = C.loaddata(catalog)       # Loading the whole catalog content.
       head = C.loadheader(catalog)     # Loading the original header.
       zps  = U.get_data(zpts,1)        # Loading Zeropoint values
       # gain = U.get_data(gains,1)     # Loading Gain Values.
       filters = get_filters(columns)   # It gets the filter names for plots.
       # Defining new variables
       ng = len(mm[:,0])                # ng is the number of galaxies.
       nf = len(mm[0,:])                # nf is the number of filters.
       if verbose: print 'ng,nl',ng,nl 
       errmag = N.zeros((ng,nf),float)  # Where the new photo errors will be saved.
       # Starting the game..
       for jj in range(nf):
           if verbose: print 'Analyzing filter %i'%(jj+1)
           # For every single band we need to read the apert_v_sigma file (area2rms)
           # to interpolate the real values of area from the catalog.
           # print 'area2rms[%i]'%(jj),area2rms
           area2rms_bands = U.get_str(area2rms,0)
           sqrtarea, sbackg, smean = U.get_data(area2rms_bands[jj],(0,1,2))
           rmsfit  = N.poly1d(N.polyfit(sqrtarea, sbackg, 3)) # CHECK
           meanfit = N.poly1d(N.polyfit(sqrtarea, smean, 2))  # CHECK
           
           # It reads the WEIGHT-map if requested
           if weight:
              if verbose: print 'Reading Weight image: ', wimas[jj] 
              wdata = fits.open(wimas[jj])[0].data
              wdatanorm = wdata / wdata.max()

           if check1:
              # Sanity plot to assure the interpolated area2sigma function was right   
              plt.figure(0, figsize = (7,6),dpi=70, facecolor='w', edgecolor='k')
              plt.clf()
              plt.plot(sqrtarea,sbackg,'ko',sqrtarea,rmsfit(sqrtarea),'r-',linewidth=2)
              plt.xlabel('$\sqrt{N}$',size=18)
              plt.ylabel('$\sigma$',size=20)
              plt.legend(['Data','Interpolation'],numpoints=1,loc='upper left') 
              plt.grid()
              plt.savefig(catalog[:-3]+'.%s.Ar2Si.check.png'%(filters[jj]),dpi=80)

              # Sanity plot to assure the interpolated area2sigma function was right   
              # plt.figure(2, figsize = (7,6),dpi=70, facecolor='w', edgecolor='k')
              plt.clf()
              plt.plot(sqrtarea, smean,'ko',sqrtarea, meanfit(sqrtarea),'r-',linewidth=2)
              plt.xlabel('$\sqrt{N}$',size=18)
              plt.ylabel('$mean$',size=20)
              plt.legend(['Data','Interpolation'],numpoints=1,loc='upper left') 
              plt.grid()
              plt.savefig(catalog[:-3]+'.%s.Ar2Mean.check.png'%(filters[jj]),dpi=80)

           fluxgal = B.mag2flux(mm[:,jj]-zps[jj]) # -meansignal+U.sqrt(fluxcorr4aperture)
           # good_sample = U.less_equal(abs(mm[:,jj]),30) # good sample
           # bad_sample  = U.greater(abs(mm[:,jj]),30)    # bad sample

           # Values to estimates the mags error.
           sqar = N.sqrt(aper)
           sigback = rmsfit(sqar)
           meansignal = meanfit(sqar)
           
           # There will be non-detected galaxies
           # with m=99 magnitudes.
           # Those numbers should not change here.
           detected     = U.less_equal(abs(mm[:,jj]),30)  # good sample
           nondetected  = U.greater(abs(mm[:,jj]),30)    # bad sample
           
           # Photom. error (as defined by SExtractor) but using the new sigma value!       
           if weight:
               pixw = N.zeros(ng)
               for hhh in range(ng): pixw[hhh]=wdatanorm[int(yy[hhh])-1,int(xx[hhh])-1]
               fluxcor = fluxgal*pixw   
               newerror=N.sqrt((aper*sigback*sigback/fluxcor**2)+(fluxcor*gain))
               newerror *= 1.0857
           else:
               newerror=N.sqrt((aper*sigback*sigback/fluxgal**2)+(fluxgal*gain))
               newerror *= 1.0857
               
           # Assesing new uncertainties.    
           errmag[detected,jj]    = newerror[detected] 
           errmag[nondetected,jj] = em[nondetected,jj]
           
           # A new figure is create to compare SExtractor vs Empirical uncert.
           if check2:       
              line = N.arange(16.,30.,0.25)
              SExline  = U.bin_stats(mm[:,jj],em[:,jj],line,stat='mean_robust') 
              aperline = U.bin_stats(mm[:,jj],errmag[:,jj],line,stat='mean_robust')   
              # plt.figure(1,figsize = (8,7),dpi=70, facecolor='w', edgecolor='k')
              plt.clf()
              plt.plot(mm[:,jj],em[:,jj],'r+',mm[:,jj],errmag[:,jj],'k+')
              plt.plot(line,SExline,'-ro',line,aperline,'-ko',linewidth=6,alpha=0.2)
              plt.legend(['$SExtractor$','$Apertures$'],numpoints=1,loc='upper left')
              plt.xlabel('$Mags$',size=17)
              plt.ylabel('$ErrMags$',size=17)
              plt.xlim(17.,30.)
              plt.ylim(0.,1.0)
              plt.grid()
              plt.savefig(catalog[:-3]+'.%s.uncert.comparison.png'%(filters[jj]),dpi=80)
              
       # The new values of mags error are now overwrited in the original data.
       vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
       data[:,evars] = errmag[:,N.arange(nf)]
       C.savedata(data,finalcat, dir="",header=head)     # Saving and creating the new catalog.
 

def create_columns_file(cluster,zs):
    """
    It creates an optimized columns file
    for the cluster.
    If zs=True it includes the Z_S parameter.
    """
    body1 = """#FILTERS  Calibration  err_zp  zp_offset  
jplus_F348.res   0,0 AB  0.1  0.00
jplus_F378.res   0,0 AB  0.1  0.00
jplus_F395.res   0,0 AB  0.1  0.00
jplus_F410.res   0,0 AB  0.1  0.00
jplus_F430.res   0,0 AB  0.1  0.00
jplus_g_sdss.res 0,0 AB  0.1  0.00
jplus_F515.res   0,0 AB  0.1  0.00
jplus_r_sdss.res 0,0 AB  0.1  0.00
jplus_F660.res   0,0 AB  0.1  0.00
jplus_i_sdss.res 0,0 AB  0.1  0.00
jplus_F861.res   0,0 AB  0.1  0.00
jplus_z_sdss.res 0,0 AB  0.1  0.00
ID		 1 	 
M_0 		 0
    """
    
    body2 = """#FILTERS  Calibration  err_zp  zp_offset  
jplus_F348.res   0,0 AB  0.1  0.00
jplus_F378.res   0,0 AB  0.1  0.00
jplus_F395.res   0,0 AB  0.1  0.00
jplus_F410.res   0,0 AB  0.1  0.00
jplus_F430.res   0,0 AB  0.1  0.00
jplus_g_sdss.res 0,0 AB  0.1  0.00
jplus_F515.res   0,0 AB  0.1  0.00
jplus_r_sdss.res 0,0 AB  0.1  0.00
jplus_F660.res   0,0 AB  0.1  0.00
jplus_i_sdss.res 0,0 AB  0.1  0.00
jplus_F861.res   0,0 AB  0.1  0.00
jplus_z_sdss.res 0,0 AB  0.1  0.00
ID		 1 	 
M_0 		 0
Z_S              0
    """
    
    if zs: body = body2
    else:  body = body1   
    filename = finalroot+'/%s/images/%s.photo.columns'%(cluster,cluster)
    fileout = open(filename,'w')
    fileout.write(body)
    fileout.close()

def replace_photo_uncert(cluster):
    """

    """
    data = C.loaddata(catalog)      # Loading the whole catalog content.
    head = C.loadheader(catalog)    # Loading the original header.
    mm = A.get_magnitudes(catalog,columns)
    em = A.get_errmagnitudes(catalog,columns)
    filters = B.get_filter_list(columns)
    
    nl = len(mm[:,0])    # nl is the number of detections inside every single band.
    nf = len(mm[0,:])    # nf is the number of bands inside the catalog. 
    errmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved. 

    for jj in range(nf):
        for ii in range(nl):
            if mm[ii,jj] == -99.: errmag[ii,jj] = 0.00
            else:  errmag[ii,jj] = em[ii,jj]   
    
    # New values of mags error overwrites now the original data.
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data[:,evars] = errmag[:,U.arange(nf)]
    C.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.


def replace_upplimits(cluster):
    """

    """
    data = C.loaddata(catalog)      # Loading the whole catalog content.
    head = C.loadheader(catalog)    # Loading the original header.
    m = get_magnitudes(catalog,columns)
    em = get_errmagnitudes(catalog,columns)
    filters = bpt.get_filter_list(columns)
    
    nl = len(m[:,0])    # nl is the number of detections inside every single band.
    nf = len(m[0,:])    # nf is the number of bands inside the catalog. 
    errmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved. 

    for jj in range(nf):
        # maglim = get_limitingmagnitude(m[:,jj],em[:,jj],3.,0.25)
        maglim = bpt.get_limitingmagnitude(m[:,jj],em[:,jj],1.,0.25)
        print 'Limiting Magnitude for filter %s: %.3f'%(filters[jj],maglim)
        for ii in range(nl):
            # print 'm[%i,%i]'%(ii,jj),m[ii,jj]
            if m[ii,jj] != 99. :         
               errmag[ii,jj] = em[ii,jj]    
            else:
               # print 'UNDETECTED OBJECT. SAVING ITS LIMITING MAGNITUDE !!'
               errmag[ii,jj] = maglim
    
    
    # New values of mags error overwrites now the original data.
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data[:,evars] = errmag[:,U.arange(nf)]
    C.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.


def get_magnitudes(catalog,columns):
    """
    It reads the magnitudes from a catalogue.
    """
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data = C.loaddata(catalog)
    mags = data[:,vars] 
    return mags


def get_errmagnitudes(catalog,columns):
    """
    It reads the magnitude uncertainties from a catalogue.
    """
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data = C.loaddata(catalog)
    emags = data[:,evars] 
    return emags


def get_usefulcolumns(columns):
 
    """
    It extracts the vars,evars,posref,zpe,zpo information
    from a columns file.
    vars & evars: mag & emags positions inside the catalog.
    ==========================================
    USAGE: 
    columns = 'cluster.columns'
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    """
    filt = U.get_str(columns,0)
    nf = 0
    for ii in range(len(filt)):
        if filt[ii][-4:] == '.res': nf += 1
        if filt[ii] == 'M_0': 
           posM0 = ii
    print 'Number of filters detected... ', nf

    filtref = int(U.get_str(cols,1)[posM0])-1
    rawvars = U.get_str(columns,1,nf)
    vars = U.zeros(nf) 
    evars = U.zeros(nf)
    for jj in range(nf):
        vars[jj] = int(rawvars[jj].split(',')[0])-1   # -1 because of columns's notation
        evars[jj] = int(rawvars[jj].split(',')[1])-1  # -1 because of columns's notation
        if vars[jj] == filtref: posref = int(vars[jj])
    zpe,zpo = U.get_data(columns,(3,4),nf)
    vars = vars.astype(int)
    evars = evars.astype(int)
    return vars,evars,posref,zpe,zpo



def match_spz_sample(cluster): # TO CHECK
       
    finalcat1 = catalog2[:-3]+'CLASH.redu.cat'
    finalcat2 = catalog2[:-3]+'nada.cat'
    # if not os.path.exists(finalcat1):
    if not os.path.exists(finalcat2):
        # print 'Final catalog does not exist yet.'                           
        if os.path.exists(catalog1) and os.path.exists(catalog2):
            # It matches up detections to its Spectroscopic Sample.
            # Reading specz catalog
            print 'Reading info1 before matching...'
            speczsample = catalog1
            idsp,xsp,ysp = U.get_data(speczsample,(0,3,4))
            goodsp = U.greater_equal(xsp,1500) * U.less_equal(xsp,3500)
            goodsp *= U.greater_equal(ysp,1500) * U.less_equal(ysp,3500)
            idsp,xsp,ysp = U.multicompress(goodsp,(idsp,xsp,ysp))
            print 'New dimension for specz catalogue: ',len(xsp)
            # rasp,decsp,xsp,ysp,zsp = get_data(speczsample,(0,1,2,3,4))
            # xsp,ysp,zsp = get_data(speczsample,(1,2,7))
            ####### idsp = U.arange(len(xsp))+1 
            # idsp = arange(len(rasp))+1
            # Reading ColorPro catalog
            print 'Reading info2 before matching...'
            idcol,xcol,ycol = U.get_data(catalog2,(0,3,4))
            print 'Dimension for input catalogue before compressing: ',len(idcol)
            gsp = U.greater_equal(xcol,1500) * U.less_equal(xcol,3500)
            gsp *= U.greater_equal(ycol,1500) * U.less_equal(ycol,3500)
            idcol,xcol,ycol = U.multicompress(gsp,(idcol,xcol,ycol))
            print 'Dimension for input catalogue after compressing: ',len(idcol)
            # Using "matching_vects" to match up samples...
            print 'Matching samples....'
            pepe = CT.matching_vects(idcol,xcol,ycol,idsp,xsp,ysp,1.1)   # We use now X,Y instead RA,Dec
            # Compressing matches for ColorPro...
            print 'Compressing matches...'
            matchidcol = pepe[:,0].astype(int)
            gdet_col = U.greater(matchidcol,0)  # Excluding 0's (non matched detections)
            matchidcol = U.compress(gdet_col,(matchidcol))
            # Compressing matches for Spectroscopic...
            matchidsp = pepe[:,1].astype(int)
            gdet_spz = U.greater(matchidsp,0)   # Excluding 0's (non matched detections)
            matchidsp = U.compress(gdet_spz,(matchidsp))
            print 'len(idcol)',len(idcol)
            print 'len(idsp)',len(idsp)
            if len(matchidcol) == len(matchidsp):
                print 'Creating idredu & zsredu '
                print 'Dimension of matchidsp ',len(matchidsp)
                idredu = U.zeros(len(matchidsp))
                idspredu = U.zeros(len(matchidsp))
                for ii in range(len(matchidsp)):
                    colindex = A.id2pos(idcol,matchidcol[ii]) # Position for Index idcol
                    spzindex = A.id2pos(idsp,matchidsp[ii])   # Position for Index idsp
                    idredu[ii] = idcol[colindex]  # ID for ColorPro
                    idspredu[ii] = idsp[spzindex]    # Specz for Specz
                    
                # A new smaller catalog will be created containing specz info as an extra column.
                print 'Selecting by rows... ' 
                finalcat1 = catalog2[:-3]+'UDF.redu.cat'
                finalcat2 = catalog2[:-3]+'CLASH.redu.cat'
                U.put_data(catalog2[:-3]+'idsfrommatch.txt',(idredu,idspredu))
                A.select_rows_bylist_sorted(catalog1,idspredu,finalcat1)
                A.select_rows_bylist_sorted(catalog2,idredu,finalcat2)               



def get_usefulcolumns(columns):
    """
    It extracts the vars,evars,posref,zpe,zpo information
    from a columns file.
    vars & evars: mag & emags positions inside the catalog. 
====USAGE=====================================
columns = 'Abell383.columns'
vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
----
    """
    filt = U.get_str(columns,0)
    nf = 0
    for ii in range(len(filt)):
        if filt[ii][-4:] == '.res': nf += 1
        if filt[ii] == 'M_0': 
           posM0 = ii

    print 'Number of filters detected... ', nf
    filtref = int(U.get_str(columns,1)[posM0])-1
    rawvars = U.get_str(columns,1,nf)
    vars = U.zeros(nf) 
    evars = U.zeros(nf)
    for jj in range(nf):
        vars[jj] = int(rawvars[jj].split(',')[0])-1   # -1 because of columns's notation
        evars[jj] = int(rawvars[jj].split(',')[1])-1  # -1 because of columns's notation
        if vars[jj] == filtref: posref = int(vars[jj])
    zpe,zpo = U.get_data(columns,(3,4),nf)
    vars = vars.astype(int)
    evars = evars.astype(int)
    return vars,evars,posref,zpe,zpo


def get_filters(columns):
    """
    It extracts the FILTERS from a columns file.
    ===========================================================
    USAGE:
columns = 'Abell383.columns'
filtros = get_filters(columns)
----
    """    
    data = U.get_str(columns,0)
    filters=[]
    for ii in range(len(data)):
        if data[ii][-4:] == '.res': 
           filters.append(data[ii])
    return filters    





def getnumfilt_observed(catalog,columns):

    """
It counts the number of filters a galaxies was not observed (m=-99)
and then it discounts that number to the total/original number of filters.
----
    """
    mags = get_magnitudes(catalog,columns)
    nf = N.shape(mags[0,:])[0]
    no = N.shape(mags[:,0])[0]
    nfobs = N.zeros(no)
    for ii in range(no):
        kk = 0
        for jj in range(nf):
               if mags[ii,jj] == -99. : kk += 1
        nfobs[ii] = nf - kk
           
    return nfobs.astype(int)    


def getnumfilt_detected(catalog,columns):

    """
It counts the number of filters a galaxies was not detected (m=+99)
and then it discounts that number to the total/original number of filters.
    """
    mags = get_magnitudes(catalog,columns)
    nf   = N.shape(mags[0,:])[0]
    no   = N.shape(mags[:,0])[0]
    nfd  = N.zeros(no)
    for ii in range(no):
        kk = 0
        for jj in range(nf):
               if mags[ii,jj] == 99. : kk += 1
        nfd[ii] = nf - kk
           
    return nfd.astype(int)    

def getJPLUSfinalids(pointing,ids):
    """
    It returns an string with
    the new ids including the pointing
    as an initial number...
    """
    idf = []
    for gg in range(len(ids)):
       ident = '%i00000'%(pointing)
       idf.append(str(int(int(ident)+ids[gg])))
    return idf


def gettingJPLUS_123limmags(clustername):
    """
    This creates a new file with an estimation
    of the 5-sigma mag-limits for 1",2",3" apertures.
    """
    sigmas = 5
    root2images  = finalroot + 'images/%s/'%(clustername)
    listimages  = root2images+'sci.list'
    imas = U.get_str(listimages,0)
    nims = len(imas)
    finalcat = root2images+'%s.photo.mlim123.cat'%(clustername)
    if not os.path.exists(finalcat):
       outfile.open(finalcat,'w')
       outfile.write('#  5-sigma threshold \n')
       outfile.write('#  mlim_1arcs  mlim_2arcs  mlim_3arcs \n')
       for ii in range(nims):
           catalog  = root2cats + get_nickname(imas[ii])+'.cat'
           f1,f2,f3,ef1,ef2,ef3,m1 = U.get_data(catalog,(23,24,25,26,27,28,21))
           s2n_1 = N.zeros(len(f1))
           s2n_2 = N.zeros(len(f2))
           s2n_3 = N.zeros(len(f3))
           
           for ggg in range(len(f1)):
               if f1[ggg]>0. and ef1[ggg]>0.: s2n_1[ggg]=f1[ggg]/(ef1[ggg]*1.)
               else: s2n_1[ggg]= -1.00  
               if f2[ggg]>0. and ef2[ggg]>0.: s2n_2[ggg]=f2[ggg]/(ef2[ggg]*1.)
               else: s2n_2[ggg]= -1.00  
               if f3[ggg]>0. and ef3[ggg]>0.: s2n_3[ggg]=f3[ggg]/(ef3[ggg]*1.)
               else: s2n_3[ggg]= -1.00  
               
           # Here we clean the sample and select the 5-sigma detections.
           good_s2n_1 = N.greater(s2n_1,sigmas-0.2) * N.less(s2n_1,sigmas+0.2)
           good_s2n_2 = N.greater(s2n_2,sigmas-0.2) * N.less(s2n_2,sigmas+0.2)
           good_s2n_3 = N.greater(s2n_3,sigmas-0.2) * N.less(s2n_3,sigmas+0.2)
           mlim1 = N.mean(m1[good_s2n_1])
           mlim2 = N.mean(m1[good_s2n_2])
           mlim3 = N.mean(m1[good_s2n_3])
           linea = '#  %.2f  %.2f  %.2f  \n'%(filtros[ii],mlim1,mlim2,mlim3)
           outfile.write(linea)


