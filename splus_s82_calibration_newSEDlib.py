__author__ = 'albertomolino'

"""
This routines is similar to "splus_s82_calibration.py" but
it uses for calibration and photo-z estimates the new library
of SEDs which has been calibrated with the S82 data without any
previous ZP-correction; i.e., using the ZPs from Laura.

It includes two runs. First using 'GOSMOSeB11recal' and then 'GOSMOSeB11_sorted_recal'
where the later is equal to the former but re-ordering the elliptical templates.

After estimating new zero-points, it will re-run BPZ on every field
so the final precision can be computed.

"""
import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U

mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/'
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root2codes = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/'

## SEDs to be used
spectra = 'COSMOSeB11new_recal'
final_root = root2cats + '%s/'%(spectra)
lista_cats = 'photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)
apertures = ['auto']
n_apers = len(apertures)
sdss_s82_spz_cat = root2cats + 'S82_SDSS_z05.cat'

n_cats = 1
########## Cross-matching with SDSS/S82.
for ggg in range(n_cats):
    splus_speczcat = cats_names[ggg][:-3] + 'spz.z05.cat'
    if not os.path.exists(splus_speczcat):
       cmd_cross_match  = "java -jar /Users/albertomolino/codigos/Stilts/stilts.jar "
       cmd_cross_match += "tmatch2 ifmt1=ascii ifmt2=ascii in1=%s "%(cats_names[ggg])
       cmd_cross_match += "in2=%s out=%s ofmt=ascii matcher=sky values1='$2 $3' "%(sdss_s82_spz_cat, splus_speczcat)
       cmd_cross_match += "values2='$1 $2' params=1 join=1and2 find=best progress=log"
       os.system(cmd_cross_match)

    ######### ZP re-calibration
    for hhh in range(n_apers):
        splus_columns = mainroot + 'splus_%s.columns'%(apertures[hhh])
        cali_cols = final_root + os.path.basename(splus_speczcat)[:-3]+'%s_cali.columns'%(apertures[hhh])

        if not os.path.exists(cali_cols):
           cmd1 = 'python %sfullcalibrator_%s.py %s '%(root2bpz,spectra,splus_speczcat)
           cmd1 += '-cols %s -outcol %s '%(splus_columns,cali_cols)
           print cmd1
           try: os.system(cmd1)
           except: continue

        # Remove temporal files created during calibration.
        cmd0 = '/bin/rm %s*temporcal*'%(root2cats)
        try: os.system(cmd0)
        except: continue

        # Save calibration plots.
        new_folder = final_root +'zp_plots/'
        if not os.path.exists(new_folder):
           cmd8 = '/bin/mkdir %s '%(new_folder)
           try: os.system(cmd8)
           except: continue

        nick_name_field = cats_names[ggg].split('/')[-1][:-4]
        new_folder_2 = new_folder +'%s/'%(nick_name_field)
        print 'new_folder_2:',new_folder_2
        if not os.path.exists(new_folder_2):
           cmd8 = '/bin/mkdir %s '%(new_folder_2)
           try: os.system(cmd8)
           except: continue

        plots_names = root2codes+'*resoffset.png'
        cmd6 = '/bin/mv %s %s'%(plots_names,new_folder_2)
        try: os.system(cmd6)
        except: continue

        if os.path.exists(cali_cols):
           bpz_cali  = final_root + os.path.basename(splus_speczcat)[:-3]+'%s_cali.bpz'%(apertures[hhh])
           flux_cali = final_root + os.path.basename(splus_speczcat)[:-3]+'%s_cali.flux_comparison'%(apertures[hhh])
           hdf5_cali = final_root + os.path.basename(splus_speczcat)[:-3]+'%s_cali.hdf5'%(apertures[hhh])

           if not os.path.exists(bpz_cali):
              cmd2  = 'python %sbpz.py %s '%(root2bpz,splus_speczcat)
              #cmd2 += '-COLUMNS %s -OUTPUT %s '%(splus_columns,bpz_cali)
              cmd2 += '-COLUMNS %s -OUTPUT %s '%(cali_cols,bpz_cali)
              cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.015 -INTERP 5 -ZMAX 0.5 '
              cmd2 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(flux_cali)
              cmd2 += '-SPECTRA %s.list -ZMIN 0.001 -PRIOR flat -DZ 0.0001 '%(spectra)
              cmd2 += '-ZP_ERRORS "0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02" '
              cmd2 += '-HDF5 yes '
              print cmd2
              os.system(cmd2)
              pausa = raw_input('paused')

        else:
            print '%s does not exist!!'%(cali_cols)
            pausa = raw_input('paused!')

        # Renaming HDF5 files
        """
        if not os.path.exists(hdf5_cali):
           temporal_hdf5 = final_root + os.path.basename(splus_speczcat).split('.')[0] + '.hdf5'
           print 'temporal_hdf5',temporal_hdf5
           print 'hdf5_cali',hdf5_cali
           cmd7 = '/bin/mv %s %s '%(temporal_hdf5,hdf5_cali)
           print cmd7
           os.system(cmd7)

        # Remove temporal files created by BPZ.
        cmd3 = '/bin/rm %s*.npy'%(root2cats)
        try: os.system(cmd3)
        except: continue

        """



"""

plt.clf()
a1,a2,a3 = plt.hist(dz[g2],basez,normed=1,alpha=0.3,facecolor='blue',linewidth=1)
a1,a2,a3 = plt.hist(dz[g5],basez,normed=1,alpha=0.3,facecolor='green',linewidth=1)
a1,a2,a3 = plt.hist(dz[g7],basez,normed=1,alpha=0.3,facecolor='red',linewidth=1)
a1,a2,a3 = plt.hist(dz[g9],basez,normed=1,alpha=0.3,facecolor='purple',linewidth=1)
plt.legend(['$Odds>0.2$ | $dz/1+z=3\%$','$Odds>0.5$ | $dz/1+z=2\%$','$Odds>0.7$ | $dz/1+z=1.5\%$','$Odds>0.9$ | $dz/1+z=1\%$'],loc='upper center',fontsize=25,numpoints=1)
a1,a2,a3 = plt.hist(dz,basez,normed=1,alpha=0.3,histtype='step',color='grey')
plt.grid()
plt.ylim(0.,50)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
plt.xlabel('$(z_{s}-z_{b})/(1+z_{s})$',size=25,labelpad=5)
plt.ylabel('$\#$',size=25,labelpad=1)

---

"""