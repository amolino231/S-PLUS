__author__ = 'albertomolino'

"""
This routine creates individual and global dz/1+z plots
from point-estimates and PDFs.

"""
import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import splus_s82_hdf5_tools as to

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
root += 'S82/Dec2017/splus_cats_NGSL/'
hdf5list = root+'hdf5.list'
bpzlist = root+'bpz/master.STRIPE82_Photometry.m21.bpz.list'
hdf5_files = U.get_str(hdf5list,0)
n_hdf5 = len(hdf5_files)
bpz_files  = U.get_str(bpzlist,0)
n_bpz = len(bpz_files)

# Sample selection.
m_max = 19.5
o_min = 0.25

for ii in range(n_bpz):
    name = os.path.basename(hdf5_files[ii])
    print name
    try: z,dp,df = to.get_PDZerrDistribution(hdf5_files[ii],bpz_files[ii],m_max,o_min)
    except: print 'Impossible to run on ',name

#to.master_PDZerrDistribution(hdf5list,bpzlist,m_max,o_min)
#print 'Impossible to run on '