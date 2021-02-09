__author__ = 'albertomolino'

#/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U

mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/'
root2newcats = root2cats + 'using_5sdss_bands/'
if not os.path.exists(root2newcats):
    cmd8 = '/bin/mkdir %s '%(root2newcats)
    os.system(cmd8)

root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
lista_cats = 'photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)
apertures = 'auto'
spectra = 'COSMOSeB11new_recal'


########## Cross-matching with SDSS/S82.
for ggg in range(n_cats):
    splus_speczcat = cats_names[ggg][:-3] + 'spz.z05.cat'
    base_name = os.path.basename(splus_speczcat)
    cali_cols = root2cats+'%s/'%(spectra)+base_name[:-3]+'%s_cali.columns'%(apertures)
    #cali_cols = splus_speczcat[:-3]+'%s_cali.columns'%(apertures)
    #print splus_speczcat
    #print cali_cols
    #pausa = raw_input('paused')

    if os.path.exists(cali_cols):
           bpz_cali  = root2newcats+base_name[:-3]+'%s_cali_5bands.bpz'%(apertures)
           flux_cali = root2newcats+base_name[:-3]+'%s_cali_5bands.flux_comparison'%(apertures)
           if not os.path.exists(bpz_cali):
              cmd2  = 'python %sbpz.py %s '%(root2bpz,splus_speczcat)
              cmd2 += '-COLUMNS %s -OUTPUT %s '%(cali_cols,bpz_cali)
              cmd2 += '-CHECK yes -SIGMA_EXPECTED 0.015 -INTERP 5 -ZMAX 0.5 '
              cmd2 += '-FLUX_COMPARISON yes -FLUX_COMPARITION %s '%(flux_cali)
              cmd2 += 'SPECTRA %s.list -ZMIN 0.001 -PRIOR flat -DZ 0.0001 '%(spectra)
              cmd2 += '-EXCLUDE "SPLUS_F0378W,SPLUS_F0395W,SPLUS_F0410W,SPLUS_F0430W,'
              cmd2 += 'SPLUS_F0515W,SPLUS_F0660W,SPLUS_F0861W" '
              cmd2 += '-ZP_ERRORS "0.03,0.02,0.02,0.02,0.02" '
              os.system(cmd2)

    else:
            print '%s does not exist!!'%(cali_cols)
            pausa = raw_input('paused!')

    # Remove temporal files created by BPZ.
    cmd3 = '/bin/rm %s*.npy'%(root2cats)
    os.system(cmd3)


