#! /usr/local/bin python    
#-*- coding: latin-1 -*-

import os,sys
import useful as U
import numpy as N
import splus_detection_threshold as sdt
from astropy.io import fits
import matplotlib.pyplot as plt

def spurious_detect_threshold(image,sexfile):
    """
    Starting with an input configuration.sex file, 
    it runs SExtractor on both sizes of the image varying 
    the threshold value across a certain range. 
    By doing this, it is possible to stablish an bottom threshold
    for which the percentage of spurious detections is not larger
    than a few percents.
======
MAKE SURE THERE IS NOT A BLANK LINE AT THE END OF THE SEx FILE !!!!

----
image = '/Volumes/amb/imagenes/f02/f02p01_F814W_4.swp.fits'
sexfile = '/Volumes/amb/catalogos/reduction_v4/f02/f02p01_colorpro_4.sex'
fspd,base = spurious_detect_threshold(image,sexfile,'yes','yes','yes','yes')

    """
    
    verbose=1

    if os.path.exists(image) and os.path.exists(sexfile):

       print 
       print 'MAKE SURE THERE IS NOT A BLANK LINE AT THE END OF THE SEx FILE !!!!'
       print
       
       min_value = 0.9
       max_value = 1.5 
       interv = 0.05 
       base = N.arange(min_value, max_value+interv, interv)
       dim = len(base)
       spd = N.zeros((dim,2),float)
       fspd = N.zeros(dim)
       newvals = N.zeros(2)
       imageinv = sdt.decapfile(image)+'_inv.fits'
 
       if not os.path.exists(imageinv):
          print 'Creating an inverse image...'
          coeff = -1.
          sdt.multiply_image_bya_number(image,coeff,imageinv)
            
       print 'Modifying SExtractor input file...'
       if os.path.exists(imageinv):
           for ii in range(dim):
               newsexcat = sdt.decapfile(sexfile)+'_thr%.2f.cat' %(base[ii])
               newsexfile = sdt.decapfile(sexfile)+'_thr%.2f.sex' %(base[ii])
               param = ['ANALYSIS_THRESH','DETECT_THRESH','CATALOG_NAME']
               newvals = [base[ii],base[ii],newsexcat]
               # newvals[0] = base[ii]
               # newvals[1] = base[ii]

               if verbose : 
                  print 'base[%i]'%(ii),base[ii]
                  print 'sexfile',sexfile
                  print 'newsexfile',newsexfile
                  print 'param',param
                  print 'newvals',newvals

               # Modifiying THRESHOLD in conf.sex.
               sdt.modifyingSExfiles(sexfile,param,newvals,newsexfile)

               print 'Running SExtractor...' 
               for ss in range(2):
                   if ss == 0: image2 = image
                   else: image2 = imageinv
                   if os.path.exists(newsexfile):
                       cmd2 =''
                       cmd2 ='sex %s -c %s' %(image2,newsexfile)
                       print cmd2
                       try: os.system(cmd2)
                       except: print 'Impossible to run SExtractor !!' 

                       print 
                       print 'Measuring detections...'
                       catout = newsexcat
                       if os.path.exists(catout):
                           # print 'YES, the catout exists...'
                           id,x,y = U.get_data(catout,(0,1,2))
                           good = N.greater(x,1500.) * N.less(x,10300.) * N.greater(y,1100.) * N.less(y,10000.)
                           id_redu = N.compress(good,id)
                           print 'Compressing the sample..'
                           spd[ii,ss] = len(id_redu)    
                       else:
                           print '%s does not exists!!'%(catout)
                           print 'Impossible to quantify percentage of spurious detections!!'
           
           print 'Estimating spurious detections...'                 
           for jj in range(dim):
               fspd[jj] = ((spd[jj,1]*1.)/(spd[jj,0]*1.))*100.   

           print 'fspd',fspd
           print 'Plotting results....'
           plt.figure(1, figsize = (12,7),dpi=80, facecolor='w', edgecolor='k')
           plt.clf()
           plt.plot(base,fspd,'-ko',linewidth=2)
           plt.plot(base,base*0.+3.,'m--',linewidth=1.5)
           plt.xlabel('Threshold ($\sigma$)'),plt.ylabel('% Spurious detections')
           plt.xlim(min_value-interv,max_value+interv)
           plt.grid()
           outname = sdt.decapfile(image)+'_thranal.png'
           plt.savefig(outname,dpi=150)

           plt.figure(2, figsize = (12,7),dpi=80, facecolor='w', edgecolor='k')
           plt.clf()
           plt.plot(base,spd[:,0],'-ko',linewidth=5)
           plt.xlabel('Threshold ($\sigma$)',size=15),plt.ylabel('Number Detected Sources',size=15)
           plt.xlim(min_value-interv,max_value+interv)
           plt.grid()
           figname = sdt.decapfile(image)+'_numdet.png'
           plt.savefig(figname,dpi=150)

           outname2 = sdt.decapfile(image)+'_thranal.txt'
           U.put_data(outname2,(fspd,base),'# fspd  base ','%.3f  %.3f')   
           plt.close()
           
           return fspd,base

    else:
        print 'Input image or input SExfile does not exist!'
        print image
        print sexfile





def decapfile(filename):
    
    """
    It purges the extension from a file
-----
file = '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg.fits'
decapfile(file) ==> '/Volumes/amb/imagenes/detections/images/calibrated/f04p01_1_acs.deg'

    """
    nick = filename.split('/')[-1:][0] 
    ending = nick.split('.')[-1:][0]
    dim = len(ending)+1
    name = filename[:-dim]

    return name

def multiply_image_bya_number(image1,number,imageout):
    data1 = fits.open(image1)[0].data
    head1 = fits.open(image1)[0].header
    fits.writeto(imageout,data1*number,head1,clobber=True)


def modifyingSExfiles(file,param,newval,outfile):

    """
    It changes an input SExtractor config. file modifying the 'param'eters
    for the inputs 'newval'ues.
    ----
    Example: 
     param = 'GAIN' // ['GAIN','FILTER','CATALOG_NAME']
     newval = 3.    // [ 3.,'Y','pepe.cat']
    """

    temp = open(file,'r')
    datos = temp.read()
    datos = datos.split('\n')
    temp.close()

    if len(param) < 5e+3: #!= len(newval): 

       try:
         dimparam = param.split(',')
         dim2 = 0
       except:
         dim2 = len(param)
 
       # print 'dim2',dim2

       dim = len(datos)
       
       v1 = ''
       v2 = ''
       
       for ii in range(dim-1): 
           t = datos[ii].split()
           t1 = t[0]
           t2 = t[1]
           # print t1,t2
           v1 += t1
           v2 += t2
           v1 += '$'
           v2 += '$'
        
       vv1 = v1.split('$')
       vv2 = v2.split('$')
       # print 'vv1', vv1[-2:]
       # print 'vv2', vv2[-2:]

       print ''
       print '==============================================' 
       if dim2 == 0:
         for ii in range(dim):
             if vv1[ii] == param:
                vv2[ii] = newval
                print 'Parameter %s updated to %s !' %(param,newval)
                # else: print 'Parameter %s did not find !!' %(param)
               
       else: 
         for hh in range(len(param)):
           # print 'len(param)',len(param)
           # print 'param',param[hh]
           # print 'newval',newval[hh]
           for ii in range(dim):
               if vv1[ii] == param[hh]: 
                  vv2[ii] = str(newval[hh])
                  print 'Parameter %s updated to %s !' %(param[hh],newval[hh])
                  # else: print 'Parameter %s did not find !!' %(param[hh])
               
       print '=============================================='
       print ''        

       output = open(outfile,'w') 
       for ss in range(dim):
           ele = '%s   %s\n' %(vv1[ss],vv2[ss])
           # print ele
           output.write(ele) 
           
       output.close()
        

    else: print 'params and newvals have different sizes!!!'
    
    



