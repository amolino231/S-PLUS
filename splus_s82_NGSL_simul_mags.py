__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import phz_plots as P
import bpz_tools as B

#roots
root_to_seds = '/Users/albertomolino/Desktop/laura/'
final_root = root_to_seds+'SimMags/'
if not os.path.exists(final_root):
    cmd8 = '/bin/mkdir %s '%(final_root)
    os.system(cmd8)
#SEDs
seds = U.get_str(root_to_seds+'seds.list',0)
nsed = len(seds)

#Filters
#filter_list='laura_splus2.list'
filter_list='laura_splus3.list'

# Redshift
z = 0. #stars

#Getting magnitudes
for ii in range(nsed):
    a,b = P.see_sed2resolution_AB(seds[ii],z,filter_list)
    aa = B.AB(a)
    newfile = final_root + 'BPZ.%s.mags.txt'%(seds[ii])
    U.put_data(newfile,(b,aa-aa[7]),'# lambda mag')


