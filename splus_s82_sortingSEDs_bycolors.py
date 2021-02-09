__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/codigos/bpz-1.99.2')
import useful as U
import numpy as N
import cosmology as C

rootSED = '/Users/albertomolino/codigos/bpz-1.99.2/SED/'
#library_name = 'NGSLnew'
library_name = 'QSO'
seds = U.get_str(rootSED+library_name+'.list',0)
nsed = len(seds)
colores = N.zeros(nsed)

#filter1 = 'HST_ACS_WFC_F435W.res'
#filter2 = 'HST_ACS_WFC_F606W.res'
filter1 = 'SPLUS_uJAVA.res'
filter2 = 'SPLUS_gSDSS.res'

for ii in range(nsed):
    sed = seds[ii]
    colores[ii] = C.color_z(sed,filter1,filter2,0.)

sortedseds,sortcols = U.multisort(colores,(seds,colores))

fileout = rootSED+library_name+'sorted.list'
newdata = open(fileout,'w')
for ii in range(nsed):
    newdata.write('%s \n'%(sortedseds[ii]))
newdata.close()
