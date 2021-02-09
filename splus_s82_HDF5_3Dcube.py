__author__ = 'albertomolino' 

import os,sys
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U
import numpy as N
import h5py
import matplotlib.pyplot as plt
import tables as tb


def collapse_HDF5_in_magnitude(hdf5file,m_min,m_max,dm):
    """
    It colapses a HDF5-file in magnitude bins.
    """
    new_hdf5_filename = hdf5file[:-4]+'compM.hdf5'
    if not os.path.exists(new_hdf5_filename):

       # Initializing the HDF5-file.
       filtros = tb.Filters(complevel=5,complib="lzo") #lz0 is much faster than zlib
       fp_file = tb.openFile(new_hdf5_filename,mode="w",title="S-PLUS HDF5 files")

       # Estimating the dimensions for the HDF5-file.
       p = h5py.File(hdf5file, mode='r')
       pdz = p.get('Likelihood')
       nz  = N.shape(pdz)[1]
       nt  = N.shape(pdz)[2]
       mm = p.get('m_0')[:]
       zz = p.get('redshift')[:]
       tt  = p.get('type')[:]
       basem = N.arange(m_min,m_max+dm,dm)
       nm = len(basem)-1  # number of magnitude-bins.

       # Including this information in the HDF5-file.
       zh = fp_file.createArray(fp_file.root,"redshift",zz)
       th = fp_file.createArray(fp_file.root,"Template",tt)
       mh = fp_file.createArray(fp_file.root,"Magnitude",basem)

       # Defining dimension for the main matrix.
       full_table = fp_file.createCArray(fp_file.root,"AllData",tb.Float32Atom(),
                    shape=(nm,nz,nt),chunkshape=(1,nz,nt),filters=filtros)

       for ss in range(nm):
           good_m_sample  = N.greater_equal(mm,basem[ss])
           good_m_sample *= N.less_equal(mm,basem[ss+1])
           #print pdz[good_m_sample,:,:]
           pdf_redu = pdz[good_m_sample,:,:]
           ng = N.shape(pdf_redu)[0]
           print '%i galaxies in %.2f<mag<%.2f '%(ng,basem[ss],basem[ss+1])
           for gg in range(ng):
               galaxy_zt_pdf = pdf_redu[gg,:,:]
               galaxy_zt_pdf /= galaxy_zt_pdf.sum()
               full_table[ss,:,:] += galaxy_zt_pdf[:,:]
           #full_table[ss,:,:] /= full_table[ss,:,:].sum()
       fp_file.close()


"""
# General roots
root = '/home/amolino/work/SPLUS/SVD/S82/Dec2018/catalogs/'
root = '/home/amolino/work/JPLUS_photoz/2018_FirstDataRelease/JPLUS_full_sample_March2018/'
# Reading HDF5 files
hdf5list = root + 'hdf5.list'
hdf5_files = U.get_str(hdf5list,0)
n_hdf5 = len(hdf5_files)

# Magnitude resolution for final file
m_min = 14
m_max = 22
dm    = 1.
base_m = N.arange(m_min,m_max+dm,dm)
n_m_ele = len(base_m)

for ss in range(n_hdf5):
    print 'Processing file: ',hdf5_files[ss]
    collapse_HDF5_in_magnitude(hdf5_files[ss],m_min,m_max,dm)    


"""
