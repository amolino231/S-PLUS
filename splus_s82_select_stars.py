__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N

#Roots
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'released_catalogues/'
final_root = root2cats + 'stars/'
if not os.path.exists(final_root):
   cmd = '/bin/mkdir %s '%(final_root)
   os.system(cmd)

# Reading list of catalogues.
list_catalogues = root2cats + 'master_SPLUS_STRIPE82_photo_BPZ.list'
cats = U.get_str(list_catalogues,0)
n_cats = len(cats)

# New file
filename = final_root + 'splus_stars.cat'
out_file = open(filename,'w')
out_file.write('#   1 ID  \n')
out_file.write('#   2 RA  \n')
out_file.write('#   3 Dec \n')
out_file.write('#   4 U_auto \n')
out_file.write('#   5 F0378_auto \n')
out_file.write('#   6 F0395_auto  \n')
out_file.write('#   7 F0410_auto \n')
out_file.write('#   8 F0430_auto  \n')
out_file.write('#   9 G_auto    \n')
out_file.write('#  10 F0515_auto  \n')
out_file.write('#  11 R_auto  \n')
out_file.write('#  12 F0660_auto  \n')
out_file.write('#  13 I_auto    \n')
out_file.write('#  14 F0861_auto  \n')
out_file.write('#  15 z_auto  \n')
out_file.write('# ID RA Dec U_auto F0378_auto F0395_auto '
               'F0410_auto F0430_auto G_auto F0515_auto '
               'R_auto F0660_auto I_auto F0861_auto z_auto'\n)

for ii in range(n_cats):
    print 'reading catalogue %i '%(ii+1)
    ids,ra,dec,u,f378,f395,f410,f430 = U.get_data(cats[ii],(0,1,2,15,24,33,42,51))
    g,f515,r,f660,i,f861,z,ps = U.get_data(cats[ii],(60,69,78,87,96,105,114,132))
    #Selecting good stars
    good  = N.greater_equal(ps,0.9)
    good *= N.less_equal(abs(u),30) * N.less_equal(abs(g),30)
    good *= N.less_equal(abs(r),30) * N.less_equal(abs(i),30)
    good *= N.less_equal(abs(z),30) * N.less_equal(abs(f378),30)
    good *= N.less_equal(abs(f395),30) * N.less_equal(abs(f410),30)
    good *= N.less_equal(abs(f430),30) * N.less_equal(abs(f515),30)
    good *= N.less_equal(abs(f660),30) * N.less_equal(abs(f861),30)
    # Compressing.
    ids,ra,dec,u,f378,f395,f410,f430 = U.multicompress(good,(ids,ra,dec,u,f378,f395,f410,f430))
    g,f515,r,f660,i,f861,z = U.multicompress(good,(g,f515,r,f660,i,f861,z))
    # Writing data into new catalogue.
    n_stars = len(ids)
    for ss in range(n_stars):
        linea = '%i  %.4f  %.4f  %.2f  '%(ids[ss],ra[ss],dec[ss],u[ss])
        linea += '%.2f  %.2f  %.2f  %.2f  '%(f378[ss],f395[ss],f410[ss],f430[ss])
        linea += '%.2f  %.2f  %.2f  %.2f  '%(g[ss],f515[ss],r[ss],f660[ss])
        linea += '%.2f  %.2f  %.2f  '%(i[ss],f861[ss],z[ss])
        out_file.write(linea+'\n')

out_file.close()