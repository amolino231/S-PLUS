__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import alhambra_photools as A
import bpz_tools as B
import splus_calib_tools as sct
import matplotlib.pyplot as plt

plots = 1

### Paths
# to data
root_to_cats = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
root_to_cats += 'S82/Dec2017/'
#root_to_flux = root_to_cats + 'splus_cats_NGSL/COSMOSeB11new_recal_OT/'

# to BPZ
root_to_bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root_to_seds = root_to_bpz+'SED/'

#final_sed_root = root_to_seds + 'newSEDcal_PriorSM/'
final_sed_root = root_to_seds + 'newSEDcal_PriorSM_May18/'
if not os.path.exists(final_sed_root):
   cmd = '/bin/mkdir %s'%(final_sed_root)
   os.system(cmd)
final_sed_root_plots = final_sed_root + 'plots/'
if not os.path.exists(final_sed_root_plots):
   cmd = '/bin/mkdir %s'%(final_sed_root_plots)
   os.system(cmd)
final_sed_root_data = final_sed_root + 'data/'
if not os.path.exists(final_sed_root_data):
   cmd = '/bin/mkdir %s'%(final_sed_root_data)
   os.system(cmd)
final_sed_root_seds = final_sed_root+ 'SEDs/'
if not os.path.exists(final_sed_root_seds):
   cmd = '/bin/mkdir %s'%(final_sed_root_seds)
   os.system(cmd)

### SED-models
#sed_library = 'COSMOSeB11new_recal.list'
sed_library = 'COSMOSeB11new.list'
sed = U.get_str(root_to_seds+sed_library,0)
n_seds = len(sed)

#Wavelength Resolution of SED corrections
new_delta_lbda = 100.
base_wavel = N.arange(2500,9000.,new_delta_lbda)
min_wave_corr = 3000
max_wave_corr = 8000
normal_wavel = 7000 # Normalization for plots

# Minimum number of galaxies to compute the corrections.
min_ng = 50

# Columns and Flux_comparison to be used.
columns = root_to_cats + 'splus_cats_NGSL/master_splus_auto.columns'
#fluxcomp  = root_to_flux + 'masterOT_COSMOSeB11new_recal_SM.flux_comparison'
fluxcomp = root_to_cats + 'splus_cats_NGSL/May18_recal_OT/masterFLUX_newZPCal.flux_comparison'

## Reading fluxes from flux_comparison
ft,fob,efob = A.get_fluxes(columns,fluxcomp)
## Reading redshift from flux_comparison
z_s = U.get_data(fluxcomp,2)
## Reading Templates from flux_comparison
t_b = U.get_data(fluxcomp,3)
## Reading magnitudes from flux_comparison
# mo= U.get_data(fluxcomp,1)
## Number of galaxies
n_gal = len(z_s)

### Getting filters and eff. wavelengths
filters = B.get_filter_list(columns)
n_filters = len(filters)
eff_waves = N.zeros(n_filters)
for ii in range(n_filters):
    eff_waves[ii] = B.effective_wavelength(filters[ii])
    print 'Filter: %s, Eff.Wav:%.3f '%(filters[ii],eff_waves[ii])

for ss in range(n_seds):
    print 'Processing SED: ',ss+1
    if ss<1:
        good_seds = N.less_equal(t_b,1.5+ss)
    elif ss==n_seds-1:
        good_seds = N.greater_equal(t_b,ss+0.5)
    else:
        good_seds = N.greater_equal(t_b,ss+0.5) * N.less_equal(t_b,ss+1.5)

    # New dimension for each SED model.
    z_s_redu = z_s[good_seds]
    new_dim = len(z_s_redu)
    print 'new_dim',new_dim
    rf_wavel = N.zeros((new_dim * n_filters),float)
    delta_f  = N.zeros((new_dim * n_filters),float)
    kk = 0
    for ii in range(n_filters):
        ft_redu = ft[ii][good_seds]
        fob_redu = fob[ii][good_seds]
        for hh in range(new_dim):
            rf_wavel[kk] = B.effective_wavelength(filters[ii])/(1.+z_s_redu[hh])
            delta_f[kk]  = ft_redu[hh]/(fob_redu[hh])
            kk += 1
            #pausa = raw_input('paused')

        clean_sample = N.less(abs(delta_f),10.0)

    average_corr = U.bin_stats(rf_wavel[clean_sample],delta_f[clean_sample],
                               base_wavel,'mean_robust')

    #This gets rid of issues with the edges
    average_corr = N.where(base_wavel<min_wave_corr,1.,average_corr)
    average_corr = N.where(base_wavel>max_wave_corr,1.,average_corr)

    if plots:
        plt.figure(1, figsize = (13.5,10.),dpi=80, facecolor='w', edgecolor='k')
        plt.clf()
        plt.subplot(211)
        plt.title('SED: %s'%(sed[ss]),size=15)
        plt.plot(rf_wavel[clean_sample],delta_f[clean_sample],'+',alpha=0.2)
        plt.plot(base_wavel,average_corr,'-ro',lw=3)
        plt.grid()
        plt.ylim(0.5 ,1.5)
        plt.xlim(min_wave_corr*0.9,max_wave_corr*1.1)
        plt.ylabel('$F_{th}/F_{ob}$',size=20)


    outfilename = final_sed_root_data + sed[ss] + 'S82SPLUS_crf_res%iAA.dat'%(new_delta_lbda)
    U.put_data(outfilename,(base_wavel,average_corr),'# base_wavel average_corr')


    ## Here it applies the corrections to the original templates
    sed_wavel, sed_flux_orig = U.get_data(root_to_seds+sed[ss],(0,1))

    corr_ori_wavel = U.match_resol(base_wavel,average_corr,sed_wavel)
    if new_dim > min_ng:
       sed_flux_new = sed_flux_orig / (1. * corr_ori_wavel)
    else:
       sed_flux_new = sed_flux_orig / 1.

    if plots:
        plt.subplot(212)
        pepe = sct.lookcloser(sed_wavel,normal_wavel)
        plt.semilogy(sed_wavel,sed_flux_orig/sed_flux_orig[pepe-1],'-',
                   linewidth=8.,alpha=0.5,color='grey')

        plt.semilogy(sed_wavel,sed_flux_new/sed_flux_new[pepe-1],'-',
                   linewidth=3.,alpha=0.7,color='red')
        plt.xlim(min_wave_corr*0.9,max_wave_corr*1.1)
        plt.ylim(0.01,10.0)
        #plt.ylim(0.0002,30.6)
        plt.grid()
        plt.xlabel('Wavelength',size=20)
        plt.ylabel('$F_{th}{new}$',size=20)

        outplot_filename = final_sed_root_plots + sed[ss] + 'S82SPLUS_crf_res%iAA.png'%(new_delta_lbda)
        plt.savefig(outplot_filename,dpi=100)

    #pausa = raw_input('paused')

    # Here it saves the new SED in the final directory
    outsed_filename = final_sed_root_seds + sed[ss][:-3] + 'recal.sed'
    U.put_data(outsed_filename,(sed_wavel,sed_flux_new),'# wavel fl')





