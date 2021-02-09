__author__ = 'albertomolino'

#/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/
"""
Here the fields are calibrated using the previous GOSMOSeB11.list library.
After ZP-recal, BPZ is run on all fields. HDF5 files are created.

"""


import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U

mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/'
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root2codes = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/'
lista_cats = 'cat/individuals/photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)
apertures = ['auto']
n_apers = len(apertures)
sdss_s82_spz_cat = root2cats + 'S82_SDSS_z1.cat'
#spectra = 'eB11'
spectra = 'COSMOSeB11new_recal'

########## Cross-matching with SDSS/S82.
for ggg in range(n_cats):
    #splus_speczcat = cats_names[ggg][:-3] + 'spz.cat'
    #splus_speczcat = cats_names[ggg][:-3] + 'spz.z02.cat'
    splus_speczcat = cats_names[ggg][:-3] + 'spz.z1.cat'
    if not os.path.exists(splus_speczcat):
       cmd_cross_match  = "java -jar /Users/albertomolino/codigos/Stilts/stilts.jar "
       cmd_cross_match += "tmatch2 ifmt1=ascii ifmt2=ascii in1=%s "%(cats_names[ggg])
       cmd_cross_match += "in2=%s out=%s ofmt=ascii matcher=sky values1='$3 $4' "%(sdss_s82_spz_cat, splus_speczcat)
       cmd_cross_match += "values2='$1 $2' params=1 join=1and2 find=best progress=log" 
       os.system(cmd_cross_match)

"""
    ######### ZP re-calibration
    for hhh in range(n_apers):
        splus_columns = mainroot + 'splus_%s.columns'%(apertures[hhh])
        #cali_cols = splus_speczcat[:-3]+'%s_%s_cali.columns'%(apertures[hhh],spectra)
        cali_cols = splus_speczcat[:-3]+'%s_cali.columns'%(apertures[hhh])
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
        new_folder = os.path.dirname(cats_names[ggg])+'/zp_plots/'
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
           bpz_cali  = splus_speczcat[:-3]+'%s_cali.bpz'%(apertures[hhh])
           flux_cali = splus_speczcat[:-3]+'%s_cali.flux_comparison'%(apertures[hhh])
           hdf5_cali = splus_speczcat[:-3]+'%s_cali.hdf5'%(apertures[hhh])
           if not os.path.exists(bpz_cali):
              cmd2  = 'python %sbpz.py %s '%(root2bpz,splus_speczcat)
              cmd2 += '-COLUMNS %s -OUTPUT %s '%(cali_cols,bpz_cali)
              #cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.015 -INTERP 5 -ZMAX 0.3 '
              cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.015 -INTERP 5 -ZMAX 0.5 '
              cmd2 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(flux_cali)
              cmd2 += '-SPECTRA %s.list -ZMIN 0.0001 -PRIOR flat -DZ 0.0001 '%(spectra)
              cmd2 += '-HDF5 yes '
              cmd2 += '-ZP_ERRORS "0.03,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02" '
              print cmd2
              os.system(cmd2)
              
        else:
            print '%s does not exist!!'%(cali_cols)
            pausa = raw_input('paused!')

        # Renaming HDF5 files


        if not os.path.exists(hdf5_cali):
           temporal_hdf5 = os.path.dirname(splus_speczcat)+'/'
           temporal_hdf5 += os.path.basename(splus_speczcat).split('.')[0] + '.hdf5'
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
field = '0027'
bpz1 = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Sept2017/splus_cats/Stripe82_%s_Photometry.spz.aper_cali.bpz'%(field)
bpz2 = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Sept2017/splus_cats/Stripe82_%s_Photometry.spz.petro_cali.bpz'%(field)
bpz3 = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Sept2017/splus_cats/Stripe82_%s_Photometry.spz.auto_cali.bpz'%(field)
bpz4 = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Sept2017/splus_cats/best.bpz';os.system('/bin/rm %s'%(bpz4))
mixing_bpz_files(bpz1,bpz2,bpz3)
z_max=0.15
a1 = B.d_stats(bpz1,zmax=z_max).nice().split('\n')[1]
a2 = B.d_stats(bpz2,zmax=z_max).nice().split('\n')[1]
a3 = B.d_stats(bpz3,zmax=z_max).nice().split('\n')[1]
a4 = B.d_stats(bpz4,zmax=z_max).nice().split('\n')[1]
print a1+'  APER'
print a2+'  PETRO'
print a3+'  AUTO'
print a4+'  MIX'

==========
field = '0021'
root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Sept2017/splus_cats/'
bpz1 = root+'eB11/Stripe82_%s_Photometry.spz.z02.auto_cali.bpz'%(field)
bpz2 = root+'GOSMOSeB11/Stripe82_%s_Photometry.spz.z02.auto_cali.bpz'%(field)
z_max=0.2;m_max=18;o_min=0.0
a1 = B.d_stats(bpz1,zmax=z_max,mmax=m_max,omin=o_min).nice().split('\n')[1]
a2 = B.d_stats(bpz2,zmax=z_max,mmax=m_max,omin=o_min).nice().split('\n')[1]
print a1+'  eB11'
print a2+'  GOSMOS'

"""