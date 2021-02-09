__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U
import numpy as N
import h5py
import matplotlib.pyplot as plt
import tables as tb

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root += 'splus_cats_NGSL/JPLUS_full_sample_March2018/hdf5/'
# Reading HDF5 files
hdf5list = root + 'hdf5.list'
hdf5_files = U.get_str(hdf5list,0)
n_hdf5 = len(hdf5_files)
minval = 1.0e-6
plots = 1


# New HDF5 filename.
new_hdf5_filename = root+'master.compM.hdf5'
if not os.path.exists(new_hdf5_filename):
   filtros = tb.Filters(complevel=5,complib="lzo")
   fp_file = tb.openFile(new_hdf5_filename,mode="w",title="Master HDF5 file")

   # Estimating the dimensions for the HDF5-file.
   p = h5py.File(hdf5_files[0], mode='r')
   pdz = p.get('AllData')
   nz  = N.shape(pdz)[1]
   nt  = N.shape(pdz)[2]
   mm = p.get('Magnitude')[:]
   nm = len(mm)
   zz = p.get('redshift')[:]
   tt  = p.get('Template')[:]

   # Including this information in the HDF5-file.
   zh = fp_file.createArray(fp_file.root,"redshift",zz)
   th = fp_file.createArray(fp_file.root,"Template",tt)
   mh = fp_file.createArray(fp_file.root,"Magnitude",mm)

   # Defining dimension for the main matrix.
   full_table = fp_file.createCArray(fp_file.root,"AllData",tb.Float32Atom(),
                shape=(nm,nz,nt),chunkshape=(1,nz,nt),filters=filtros)


for ii in range(n_hdf5):
    try: 
         p = h5py.File(hdf5_files[ii], mode='r')
         pepe = p.get('AllData')
         print 'reading file %i/%i: %s'%(ii+1,n_hdf5,os.path.basename(hdf5_files[ii]))
         if ii<1:
            zz = p.get('redshift')[:]
            nz = len(zz)
            tt = p.get('Template')[:]
            nt = len(tt)
            m = p.get('Magnitude')[:]
            nm = len(m)
            final_mat = N.zeros((nm,nz,nt),float)

         for ss in range(nm-1):
             if ii<1:
                 full_table[ss,:,:] = N.where(pepe[ss,:,:]>minval,pepe[ss,:,:]/(pepe[ss,:,:].max()),0)
                 final_mat[ss,:,:]=N.where(pepe[ss,:,:]>minval,pepe[ss,:,:]/(pepe[ss,:,:].max()),0)
             else:
                 full_table[ss,:,:] += N.where(pepe[ss,:,:]>minval,pepe[ss,:,:]/(pepe[ss,:,:].max()),0)
                 final_mat[ss,:,:]+=N.where(pepe[ss,:,:]>minval,pepe[ss,:,:]/(pepe[ss,:,:].max()),0)


    except:
         print 'Impossible to read: ', os.path.basename(hdf5_files[ii])



#for ss in range(nm):
#    full_table[ss,:,:] = final_mat[ss,:,:]
fp_file.close()



#if plots:
plt.figure(22,figsize = (13,6),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
for ss in range(4):
       plt.subplot(1,4,ss+1)
       if ss<1:
          yo = U.sum(final_mat[0:3,50:-100,:],axis=0)
       else:
          yo = final_mat[3+ss,50:-100,:]
       yo2 = N.where(yo<minval,0.,yo)
       print yo2.min()
       print yo2.max()
       print ''
       plt.contour(zz[50:-100],tt,N.log10(yo2).T,800,linewidths=2,
                   vmin=-1.0,vmax=1.)
       if ss<1:
           plt.title('R$\leq$%s'%(m[ss+3]),size=20)
       else:
           plt.title('R=%s'%(m[ss+3]),size=20)
       if ss in [0,4]:
          plt.ylabel('Spectral-type',size=20)
       #if ss in [4,5,6,7]:
       plt.xlabel('$z$',size=25,labelpad=-2)
plt.savefig('/Users/albertomolino/Desktop/hdf5plot2.png',dpi=90)

basez = N.arange(0.005,0.45,0.03)
plt.figure(23,figsize = (13,6),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
for ss in range(4):
       plt.subplot(1,4,ss+1)
       pepe0=final_mat[ss,50:-150,:]
       pepa0=U.sum(pepe0,axis=1)
       new_pdz0=U.match_resol(zz[50:-150],pepa0,basez)
       plt.plot(basez,(new_pdz0/new_pdz0.sum()*1.),'-ro')
       if ss<1:
           plt.title('R$\leq$%s'%(m[ss+3]),size=20)
       else:
           plt.title('R=%s'%(m[ss+3]),size=20)
       if ss in [0,4]:
          plt.ylabel('Spectral-type',size=20)
       plt.xlabel('$z$',size=25,labelpad=-2)
       plt.xlim(0.005,0.35)
plt.savefig('/Users/albertomolino/Desktop/hdf5plot3.png',dpi=90)

"""
ii=0



"""




"""
   plt.clf()
   for ss in range(nm-1):
       plt.subplot(2,4,ss+1)
       #plt.figure(ss+1)
       yo = full_table[ss,:,:]
       yo2 = N.where(yo<minval,0.,yo)
       plt.contour(zz,tt,N.log10(yo2).T,800,linewidths=2)
       plt.title('%s'%(m[ss]),size=20)
       if ss in [0,4]:
          plt.ylabel('Spectral-type',size=20)
       if ss in [4,5,6,7]:
          plt.xlabel('$z$',size=30)
   plt.savefig('/Users/albertomolino/Desktop/hdf5plot.png',dpi=90)


"""


"""
for ss in range(7):
     plt.figure(ss+1)
     #plt.subplot(2,4,ss+1)
     for ii in range(n_hdf5):
         p = h5py.File(lista[ii], mode='r')
         pepe = p.get('AllData')
         if ii<1:
             yo = pepe[ss,:,:]/(pepe[ss,:,:].max())
         else:    
             yo += pepe[ss,:,:]/(pepe[ss,:,:].max())
     yo2 = N.where(yo<minval,0.,yo)        
     plt.contour(zz,tt,N.log10(yo2).T,800,linewidths=2)
     plt.title('%s'%(m[ss]),size=20)

"""


