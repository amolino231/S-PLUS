__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N

"""

"""


mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/'
lista_cats = 'photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)
apertures = ['auto']
n_apers = len(apertures)
root2bpzs = root2cats + 'COSMOSeB11new_recal_OT_after_SEDrecal/'

final_root_odds = root2bpzs+'Odds/'
if not os.path.exists(final_root_odds):
    cmd8 = '/bin/mkdir %s '%(final_root_odds)
    os.system(cmd8)

new_filename = final_root_odds + 'Odds_quality.txt'
filename = open(new_filename,'w')
filename.write('# Field NumGal  <Odds> \n')


for ii in range(n_cats):
    splus_speczcat = cats_names[ii][:-3] + 'spz.z05.cat'
    nickname = os.path.basename(splus_speczcat)
    bpz_cali  = root2bpzs + nickname[:-3]+'%s_cali.OT.bpz'%(apertures[0])
    zb,mo,od,zs,chi2,tb = U.get_data(bpz_cali,(1,10,5,9,8,4))
    #x,y,r= U.get_data(cats_names[ii],(3,4,78))
    #red = N.less_equal(tb,2.5) * N.greater(zs,0.0) * N.less_equal(zs,0.3) * N.less(abs(r),19.)
    good_z = N.less_equal(zs,0.5) * N.less(abs(mo),20.)
    filename.write('%s  %i  %.3f \n'%(nickname,len(od[good_z]),N.mean(od[good_z])))
    print '%s, #: %i, <Odds>: %.3f '%(nickname,len(od[good_z]),N.mean(od[good_z]))
filename.close()

