__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import splus_calib_tools as sct
import matplotlib.pyplot as plt

### It would be interesting to extract from the images
### the exposure times per filter and include it in the figure.

# filter names
filters = ['U','F378','F395','F410','F430','G',
          'F515','R','F660','I','F861','Z']

#colors for different filters
colores = N.zeros((3,12),float)
colores[:,0]=(0.00,0.00,1.00)
colores[:,1]=(0.00,0.25,1.00)
colores[:,2]=(0.00,0.65,1.00)
colores[:,3]=(0.00,0.50,0.00)
colores[:,4]=(0.85,0.65,0.00)
colores[:,5]=(0.75,0.50,0.00)
colores[:,6]=(0.80,0.25,0.00)
colores[:,7]=(1.00,0.00,0.00)
colores[:,8]=(0.85,0.00,0.00)
colores[:,9]=(0.65,0.00,0.00)
colores[:,10]=(0.35,0.00,0.00)
colores[:,11]=(0.25,0.00,0.00)


root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root_to_cats = root+'splus_cats/'
root_to_ind  = root+'splus_indiv_cats/'
root_to_zps  = root_to_cats + 'ZPs/'
final_path = root + 'data_quality/'
lista_cats = root_to_cats+'photometry.list'
cats = U.get_str(lista_cats,0)
n_cats = len(cats)
n_fields = 4 # variables to plot: FWHM, b/a, ZPs,

# Reading exptime and backg.
exposures = final_path + 'info_exptime.txt'
backg_matrix = final_path + 'info_backg.txt'
# Reading exposure-time and background from external files.
backg_1pix = U.get_2Darray(backg_matrix)
exposuretime = U.get_2Darray(exposures)

#Position in the catalogue for each filter
mag_pos = N.array([ 18,  27,  36,  45,  54,  63,  72,  81,  90,  99, 108, 117])
s2n_pos = N.array([ 20,  29,  38,  47,  56,  65,  74,  83,  92,  101, 110, 119])

#basis
base_fwhm  = N.arange(0.5,3.5,0.15)
base_zps   = N.arange(18,23,0.25)
base_ba    = N.arange(0.8,1.01,0.01)
base_backg = N.arange(0.0,0.2,0.01)*10.
base_exp   = N.arange(0,20.,1.)
base_mlim  = N.arange(15,30,0.25)
base_mlim2 = base_mlim[:-1]+((base_mlim[1]-base_mlim[0])/2.)

## Defining variables where to store info.
fwhm_values   = N.zeros((n_cats,12),float)
boa_values    = N.zeros((n_cats,12),float)
zps_values    = N.zeros((n_cats,12),float)
mlim_values   = N.zeros((n_cats,12),float)
#seeing_values = N.zeros((n_cats,12),float)

## Reading data from catalogues
for sss in range(n_cats):
    field = os.path.basename(cats[sss])[9:-15]
    master_cat = root_to_cats + 'STRIPE82-%s_Photometry.cat'%(field) ##
    print 'reading catalog %i out of %i '%(sss+1,n_cats)
    fwhm,mr = U.get_data(master_cat,(8,84))
    fwhm_master,stars = sct.get_seeing_from_data_pro(fwhm,mr)
    for hhh in range(12):
        ind_catalog = root_to_ind+'sex_STRIPE82-%s_%s_swp.cat'%(field,filters[hhh]) ##
        fwhm_ind  = U.get_data(ind_catalog,6)
        fwhm_redu = fwhm_ind[stars]
        good_stars = N.greater_equal(fwhm_redu,1.)
        fwhm_values[sss,hhh] = U.mean_robust(fwhm_redu[good_stars])*0.55
        # aa & bb have to be picked from other catalogues.
        #aa,bb = U.get_data(ind_catalog,(11,12))
        #pausa = raw_input('paused in filter: %s'%(filters[hhh]))
        #boa_values[sss,hhh] = U.mean_robust((bb/aa)[stars])

        # Reading s/n from master for each filter
        mp,s2n = U.get_data(master_cat,(mag_pos[hhh],s2n_pos[hhh]))
        mp += 0.083
        c0 = N.greater_equal(s2n,3.) * N.less(mp,30)
        w1,w2 = N.histogram(mp[c0],base_mlim)
        mlim_values[sss,hhh] = base_mlim2[N.argmax(w1)]

    # Reading ZP values
    zp_file = root_to_zps+'STRIPE82-%s_ZP.cat'%(field)
    zps_values[sss,:] = U.get_data(zp_file,1)[:]



# Super PLOT
plt.figure(10,figsize = (12,10),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
# First we fill in with FWHM values
for iii in range(12):
    # FWHM
    plt.subplot(12,n_fields,(n_fields*iii)+1)
    fff = fwhm_values[:,iii]
    fff2 = fff[fff>0.]
    w1,w2,w3 = plt.hist(fff2,base_fwhm,
                        color= colores[:,iii],alpha=0.5,normed=1)

    plt.grid()
    plt.xlim(0.51,2.99)
    plt.ylabel('%s'%(filters[iii]),size=14.5)
    plt.yticks([])
    if iii == 11: plt.xlabel('$seeing$ $[arcs]$',size=20,labelpad=15)

    # Backg
    plt.subplot(12,n_fields,(n_fields*iii)+2)
    w1,w2,w3 = plt.hist(backg_1pix[:,iii]*10.,base_backg,
                        color= colores[:,iii],alpha=0.5,normed=1)
    plt.xlim(0.01,1.9)
    plt.grid()
    plt.yticks([])
    if iii == 11: plt.xlabel('$\sigma_{backg}^{1pix}$ $[10^{-2}]$',size=20,labelpad=15)

    # EXPTIME
    plt.subplot(12,n_fields,(n_fields*iii)+3)
    w1,w2,w3 = plt.hist(exposuretime[:,iii]/60.,base_exp,
                        color= colores[:,iii],alpha=0.5,normed=1)
    plt.xlim(.1,15.9)
    plt.grid()
    plt.yticks([])
    if iii == 11: plt.xlabel('$T_{exp}$ $[min]$',size=20,labelpad=15)

    # Mlim
    plt.subplot(12,n_fields,(n_fields*iii)+4)
    w1,w2,w3 = plt.hist(mlim_values[:,iii],base_mlim,
                        color= colores[:,iii],alpha=0.5,normed=1)
    plt.xlim(18.1,22.9)
    plt.grid()
    plt.yticks([])
    if iii == 11: plt.xlabel('$m_{comp}^{s/n>3}$',size=20,labelpad=15)


