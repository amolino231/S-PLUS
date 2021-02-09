__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U

"""
This routine is executed to run BPZ on its OT=yes mode
and obtained the expected vs observed fluxes for all galaxies
with spec-z info. These files will be used then to compute re-calibrations
to the initial SED library using 'splus_SEDrecal.py'.
Finally, new ZP-corrections will be performed for the new SED library
and BPZ will run again over all catalogues.

"""

mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/'
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root2codes = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/'
lista_cats = 'cat/individuals/photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)
lista_columns = root2cats + 'PriorSM/masterCOLUMNS_PriorSM.list'
cols_names = U.get_str(lista_columns,0)
n_cols = len(cols_names)

if n_cats <> n_cols:
    print 'Dimensions mismatch'
    sys.exit()

apertures = ['auto']
n_apers = len(apertures)

### spectral library
spectra = 'COSMOSeB11new_recal'

final_root = root2cats+'%s_OT/'%(spectra) # New root
#final_root = root2cats+'%s_OT_after_SEDrecal/'%(spectra) # New root
if not os.path.exists(final_root):
    cmd8 = '/bin/mkdir %s '%(final_root)
    os.system(cmd8)

# splus_columns = root2cats + 'splus11_auto.columns'
#sdss_s82_spz_cat = root2cats + 'S82_SDSS_z05.cat'
sdss_s82_spz_cat = root2cats + 'S82_SDSS_z1.cat'

########## Cross-matching with SDSS/S82.
for ggg in range(n_cats):
    splus_speczcat = cats_names[ggg][:-3] + 'spz.z1.cat'
    #splus_speczcat = cats_names[ggg][:-3] + 'spz.z05.cat'
    for hhh in range(n_apers):
        splus_columns = root2cats + 'master_splus_auto.columns'
        nickname = os.path.basename(splus_speczcat)
        cali_cols = cols_names[ggg]
        #cali_cols = root2cats + 'columns/'
        #cali_cols += nickname[:-3]+'%s_cali.columns'%(apertures[hhh])
        if os.path.exists(cali_cols):
           bpz_cali  = final_root+nickname[:-3]+'%s_cali.OT.bpz'%(apertures[hhh])
           flux_cali = final_root+nickname[:-3]+'%s_cali.OT.flux_comparison'%(apertures[hhh])
           #hdf5_cali = final_root+nickname[:-3]+'%s_cali.OT.hdf5'%(apertures[hhh])
           if not os.path.exists(bpz_cali):
              cmd2  = 'python %sbpz.py %s '%(root2bpz,splus_speczcat)
              cmd2 += '-COLUMNS %s -OUTPUT %s '%(splus_columns,bpz_cali)

              ### Here i am using the original zero-points.

              #cmd2 += '-COLUMNS %s -OUTPUT %s '%(cali_cols,bpz_cali)
              cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.015 -INTERP 5 -ZMAX 1.0 '
              cmd2 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(flux_cali)
              cmd2 += '-SPECTRA %s.list -ZMIN 0.005 -PRIOR SM -DZ 0.001 '%(spectra)
              cmd2 += '-HDF5 no -ONLY_TYPE yes '
              cmd2 += '-ZP_ERRORS "0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02" '
              os.system(cmd2)


        """
        # Renaming HDF5 files
        if not os.path.exists(hdf5_cali):
           temporal_hdf5 = os.path.dirname(splus_speczcat)+'/'
           temporal_hdf5 += os.path.basename(splus_speczcat).split('.')[0] + '.hdf5'
           print 'temporal_hdf5',temporal_hdf5
           print 'hdf5_cali',hdf5_cali
           cmd7 = '/bin/mv %s %s '%(temporal_hdf5,hdf5_cali)
           print cmd7
           os.system(cmd7)
        """

        # Remove temporal files created by BPZ.
        cmd3 = '/bin/rm %s*.npy'%(root2cats)
        try: os.system(cmd3)
        except: continue





"""
# For every field, i want to compare the precision
# reached by any of the three apertures.
for ggg in range(n_cats):
    print 'Analyzing field: ',os.path.basename(cats_names[ggg])
    print '================================================================================'
    print '  med    std_mad  std_ci std_phat  std  n>5sigma  n>5.*0.0300 num. '
    for hhh in range(n_apers):
        bpz_cali  = cats_names[ggg][:-3] + 'spz.%s_cali.bpz'%(apertures[hhh])
        if os.path.exists(bpz_cali):
           try:
                z_max = 8.2
                m_max = 47
                o_min = 0.9
                n_gal = len(U.get_data(bpz_cali,0))
                aa = B.d_stats(bpz_cali,zmax=z_max,mmax=m_max,omin=o_min).nice().split('\n')[1]
                bb = B.d_stats(bpz_cali,zmax=z_max,mmax=m_max,omin=0.5).nice().split('\n')[1]
                cc = B.d_stats(bpz_cali,zmax=z_max,mmax=m_max,omin=0.0).nice().split('\n')[1]
                n_gal_odds = float(bb.split()[-1])
                print aa,apertures[hhh],'o>%.2f'%(o_min)
                print bb,apertures[hhh],'o>0.5','%i'%((100.*n_gal_odds)/n_gal)
                print cc,apertures[hhh],'o>0.'
                print ' '
           except:
                print 'Not enough sample. Skipping this field.'
    print '================================================================================'
    print ' '
    pausa = raw_input('paused')
"""

"""
import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import bpz_tools as B
import useful as U
bpzlist = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
bpzlist += 'splus_cats_NGSL/bpz/master.STRIPE82_Photometry.m21.bpz.list'
bpz_files  = U.get_str(bpzlist,0)
n_bpzs = len(bpz_files)
med_value = N.zeros((n_bpzs,18),float)
std_value = N.zeros((n_bpzs,18),float)
out_value = N.zeros((n_bpzs,18),float)
num_sourc = N.zeros((n_bpzs,18),float)

for ii in range(n_bpzs):
    print bpz_files[ii]
    ao = B.d_stats(bpz_files[ii],mmin=16,mmax=18,zmax=0.5)
    pepe = ao.types().split('\n')[1:-1]
    for jj in range(18):
        pepa = pepe[jj].split()
        #print pepa
        med_value[ii,jj] = float(pepa[1])
        std_value[ii,jj] = float(pepa[2])
        out_value[ii,jj] = float(pepa[3])
        num_sourc[ii,jj] = float(pepa[4])

print 'Model  med  std  out  num '
for jj in range(18):
    linea  = '%i, %.3f, %.3f, '%(jj+1, U.mean_robust(med_value[:,jj]),U.mean_robust(std_value[:,jj]))
    linea += '%.3f  %i '%(U.mean_robust(out_value[:,jj]),U.sum(num_sourc[:,jj]))
    print linea

=== results ===
Model  med  std  out  num
1,   0.010, 0.020, 0.000  420 *
2,   0.012, 0.016, 0.000  526 *
3,   0.018, 0.021, 0.000  781 *
4,   0.019, 0.030, 0.000  706 *
5,  -0.007, 0.042, 0.000  1191 *
6,  -0.003, 0.030, 0.000  1042 *
7,  -0.002, 0.030, 0.000  1275 *
8,  -0.002, 0.026, 0.000  1499
9,  -0.003, 0.021, 0.000  1989
10, -0.003, 0.016, 0.000  2569
11, -0.003, 0.018, 0.000  810
12, -0.008, 0.027, 0.000  481
13, -0.005, 0.021, 0.000  351
14, -0.017, 0.030, 0.000  376 *
15,  0.003, 0.036, 0.037  2631 *
16,  0.004, 0.029, 0.019  3501 *
17,  0.010, 0.033, 0.000  1672 *
18,  0.002, 0.006, 0.000  282

base_sz = N.arange(0.0005,0.1,0.015)
base_sz_2 = N.arange(-0.1,0.1,0.015)
plt.clf()
for ss in range(18):
    plt.subplot(3,6,ss+1)
    good = N.greater(std_value[:,ss],0.001)
    a1,a2,a3 = plt.hist(std_value[good,ss],base_sz,alpha=0.5,normed=1)
    b1,b2,b3 = plt.hist(med_value[:,ss],base_sz_2,facecolor='red',alpha=0.5,normed=1)
    plt.legend(['S: %.3f'%(U.mean(std_value[:,ss])),'M:%.3f'%(U.mean_robust(med_value[good,ss]))]
    ,loc='upper left',fontsize=10)
    plt.grid()

"""