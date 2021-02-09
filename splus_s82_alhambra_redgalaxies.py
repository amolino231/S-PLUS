__author__ = 'albertomolino'

#/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/

# First of all, we need to recalibrate the data according to our new library of models
# and compare the results we get for those red galaxies with the originals from ALH.

# IDEALLY, it would be necessary to calibrate CCD by CCD with
# the new library so we make a fair comparison with the GOLD-catalogue.

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')

# Sample of galaxies to be used.
#sample = 'red'
sample = 'all'

# roots
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
mainroot += 'splus_cats_NGSL/ALHAMBRA_%s_galaxies/'%(sample)
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'

# BPZ library
#spectra = 'COSMOSeB11new_recal'
spectra = 'eB11'

# catalogues.
alh_cat     = mainroot + 'alhambra_%s_galaxies_specz.cat'%(sample)
alh_columns = mainroot + 'alhambra_%s_galaxies_specz.columns'%(sample)
alh_calcols = mainroot + 'alhambra_%s_galaxies_specz_recal_%s.columns'%(sample,spectra)
alh_bpz     = mainroot + 'alhambra_%s_galaxies_specz_recal_%s.bpz'%(sample,spectra)
alh_flux    = mainroot + 'alhambra_%s_galaxies_specz_recal_%s.flux_comparison'%(sample,spectra)

if not os.path.exists(alh_calcols):
    cmd1 = 'python %sfullcalibrator_%s_alhambra.py %s '%(root2bpz,spectra,alh_cat)
    cmd1 += '-cols %s -outcol %s '%(alh_columns,alh_calcols)
    print cmd1
    os.system(cmd1)

if not os.path.exists(alh_bpz):
   cmd2  = 'python %sbpz.py %s '%(root2bpz,alh_cat)
   cmd2 += '-COLUMNS %s -OUTPUT %s '%(alh_calcols,alh_bpz)
   cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.010 -INTERP 5 -ZMAX 7.0 '
   cmd2 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(alh_flux)
   cmd2 += 'SPECTRA %s.list -ZMIN 0.001 -PRIOR flat -DZ 0.0001 '%(spectra)
   cmd2 += '-ZP_ERRORS "0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,'
   cmd2 += '0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03,0.03" '
   os.system(cmd2)

"""
sample = 'red'
ruta = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/'
ruta += 'Dec2017/splus_cats_NGSL/ALHAMBRA_%s_galaxies/'%(sample)
b0 = ruta + 'alhambra_%s_galaxies_specz_recal_eB11.bpz'%(sample)
b1 = ruta + 'alhambra_%s_galaxies_specz_recal_COSMOSeB11new_recal.bpz'%(sample)
ruta2 = '/Users/albertomolino/doctorado/photo/catalogos/specz/'
b2 = ruta2+'ALHAMBRA.spz.cat'
zb0,zs0,m0,t0 = U.get_data(b0,(1,9,10,4))
dz0 = (zb0-zs0)/(1.+zs0)
zb1,zs1,m1,t1 = U.get_data(b1,(1,9,10,4))
dz1 = (zb1-zs1)/(1.+zs1)
#zb2,zs2,m2,t2 = U.get_data(b2,(72,80,65,75))
zb2,zs2,m2,t2 = U.get_data(b2,(1,9,10,4))
dz2 = (zb2-zs2)/(1.+zs2)
basem = N.arange(14.5,19.5,0.5)
values = N.zeros(len(basem))
values0 = N.zeros(len(basem))
values1 = N.zeros(len(basem))
values2 = N.zeros(len(basem))
for ii in range(len(basem)):
    g0=N.less(m0,basem[ii])
    g1=N.less(m1,basem[ii])
    #g2=N.less(m2,basem[ii])
    values[ii]=U.std_mad(dz0[g0])/(1.*U.std_mad(dz1[g1]))
    values0[ii]=U.std_mad(dz0[g0])
    values1[ii]=U.std_mad(dz1[g1])
    #values2[ii]=U.std_mad(dz2[g2])
plt.clf()
plt.subplot(211)
plt.plot(basem,values,'-ko')
plt.subplot(212)
plt.plot(basem,values0,'b--',basem,values1,'r--',lw=6)
#plt.plot(basem,values0,'b--',basem,values1,'r--',basem,values2,'k--',lw=6)

"""