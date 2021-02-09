__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N

# Paths and S-PLUS catalogues.
root2cats = '/Users/albertomolino/Desktop/SPLUS_cats/ZPs_SPLUS_to_SDSS_S82/'
lista_cats = 'photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)

# SDSS_S82_star_catalogue to compare with.
sdss_s82_spz_cat = root2cats + 'SDSS_S82_star.cat'

########## Cross-matching with SDSS/S82.
for ggg in range(n_cats):    
    splus2sdss_cat =  root2cats+os.path.basename(cats_names[ggg])[:-4]+'_to_sdsss82.cat'
    if not os.path.exists(splus2sdss_cat):
       cmd_cross_match  = "java -jar /Users/albertomolino/codigos/Stilts/stilts.jar "
       cmd_cross_match += "tmatch2 ifmt1=ascii ifmt2=ascii in1=%s "%(cats_names[ggg])
       cmd_cross_match += "in2=%s out=%s ofmt=ascii matcher=sky values1='$2 $3' "%(sdss_s82_spz_cat, splus2sdss_cat)
       cmd_cross_match += "values2='$2 $3' params=1 join=1and2 find=best progress=log"
       print cmd_cross_match
       os.system(cmd_cross_match)
       
    ### Creating an Output file where to store ZPs.
    out_filename = splus2sdss_cat[:-3]+'ZPs.info.txt'
    filename = open(out_filename,'w')
    # Header
    filename.write('# FILTER PETRO AUTO CIRC \n')
       
    ## Reading data from catalogues (3 apertures! + 1 SDSS)
    ## PETRO
    u_petro,g_petro,r_petro,i_petro,z_petro = U.get_data(splus2sdss_cat,(18,63,81,99,117))
    ## AUTO
    u_auto,g_auto,r_auto,i_auto,z_auto = U.get_data(splus2sdss_cat,(15,60,78,96,114))    
    ## CIRCULAR-3"
    u_aper,g_aper,r_aper,i_aper,z_aper = U.get_data(splus2sdss_cat,(21,66,84,102,120))
    ## SDSS
    u_sdss,g_sdss,r_sdss,i_sdss,z_sdss = U.get_data(splus2sdss_cat,(126,127,128,129,130))
    
    ## Cleaning samples (I)
    good_u_petro = N.greater(u_petro,0) * N.less(u_petro,30) 
    good_g_petro = N.greater(g_petro,0) * N.less(g_petro,30)
    good_r_petro = N.greater(r_petro,0) * N.less(r_petro,30)
    good_i_petro = N.greater(i_petro,0) * N.less(i_petro,30)
    good_z_petro = N.greater(z_petro,0) * N.less(z_petro,30)
    #
    good_u_auto = N.greater(u_auto,0) * N.less(u_auto,30) 
    good_g_auto = N.greater(g_auto,0) * N.less(g_auto,30)
    good_r_auto = N.greater(r_auto,0) * N.less(r_auto,30)
    good_i_auto = N.greater(i_auto,0) * N.less(i_auto,30)
    good_z_auto = N.greater(z_auto,0) * N.less(z_auto,30)
    #
    good_u_aper = N.greater(u_aper,0) * N.less(u_aper,30) 
    good_g_aper = N.greater(g_aper,0) * N.less(g_aper,30)
    good_r_aper = N.greater(r_aper,0) * N.less(r_aper,30)
    good_i_aper = N.greater(i_aper,0) * N.less(i_aper,30)
    good_z_aper = N.greater(z_aper,0) * N.less(z_aper,30)    
    #
    good_u_sdss = N.greater(u_sdss,0) * N.less(u_sdss,30) 
    good_g_sdss = N.greater(g_sdss,0) * N.less(g_sdss,30)
    good_r_sdss = N.greater(r_sdss,0) * N.less(r_sdss,30)
    good_i_sdss = N.greater(i_sdss,0) * N.less(i_sdss,30)
    good_z_sdss = N.greater(z_sdss,0) * N.less(z_sdss,30) 
    
    ### Computing offsets
    ### U-band
    dm_u_auto = u_auto[good_u_sdss*good_u_auto]-u_sdss[good_u_sdss*good_u_auto]
    offset_dm_u_auto = U.mean_robust(dm_u_auto) 
    dm_u_petro = u_petro[good_u_sdss*good_u_petro]-u_sdss[good_u_sdss*good_u_petro]
    offset_dm_u_petro = U.mean_robust(dm_u_petro)
    dm_u_aper = u_auto[good_u_sdss*good_u_aper]-u_sdss[good_u_sdss*good_u_aper]
    offset_dm_u_aper = U.mean_robust(dm_u_aper) 
    # Saving info.
    line  = 'uJAVA  %.3f  %.3f '%(offset_dm_u_auto,offset_dm_u_petro)
    line +=' %.3f \n'%(offset_dm_u_aper)
    filename.write(line)
    
    ### G-band
    dm_g_auto = g_auto[good_g_sdss*good_g_auto]-g_sdss[good_g_sdss*good_g_auto]
    offset_dm_g_auto = U.mean_robust(dm_g_auto)  
    dm_g_petro = g_petro[good_g_sdss*good_g_petro]-g_sdss[good_g_sdss*good_g_petro]
    offset_dm_g_petro = U.mean_robust(dm_g_petro)
    dm_g_aper = g_auto[good_g_sdss*good_g_aper]-g_sdss[good_g_sdss*good_g_aper]
    offset_dm_g_aper = U.mean_robust(dm_g_aper) 
    # Saving info.
    line  = 'gSDSS  %.3f  %.3f '%(offset_dm_g_auto,offset_dm_g_petro)
    line +=' %.3f \n'%(offset_dm_g_aper)
    filename.write(line)
    
    ### R-band
    dm_r_auto = r_auto[good_r_sdss*good_r_auto]-r_sdss[good_r_sdss*good_r_auto]
    offset_dm_r_auto = U.mean_robust(dm_r_auto) 
    dm_r_petro = r_petro[good_r_sdss*good_r_petro]-r_sdss[good_r_sdss*good_r_petro]
    offset_dm_r_petro = U.mean_robust(dm_r_petro)
    dm_r_aper = r_auto[good_r_sdss*good_r_aper]-r_sdss[good_r_sdss*good_r_aper]
    offset_dm_r_aper = U.mean_robust(dm_r_aper) 
    # Saving info.
    line  = 'rSDSS  %.3f  %.3f '%(offset_dm_r_auto,offset_dm_r_petro)
    line +=' %.3f \n'%(offset_dm_r_aper)
    filename.write(line)
    
    ### I-band
    dm_i_auto = i_auto[good_i_sdss*good_i_auto]-i_sdss[good_i_sdss*good_i_auto]
    offset_dm_i_auto = U.mean_robust(dm_i_auto) 
    dm_i_petro = i_petro[good_i_sdss*good_i_petro]-i_sdss[good_i_sdss*good_i_petro]
    offset_dm_i_petro = U.mean_robust(dm_i_petro)
    dm_i_aper = i_auto[good_i_sdss*good_i_aper]-i_sdss[good_i_sdss*good_i_aper]
    offset_dm_i_aper = U.mean_robust(dm_i_aper) 
    # Saving info.
    line  = 'iSDSS  %.3f  %.3f '%(offset_dm_i_auto,offset_dm_i_petro)
    line +=' %.3f \n'%(offset_dm_i_aper)
    filename.write(line)
    
    ### Z-band
    dm_z_auto = z_auto[good_z_sdss*good_z_auto]-z_sdss[good_z_sdss*good_z_auto]
    offset_dm_z_auto = U.mean_robust(dm_z_auto) 
    dm_z_petro = z_petro[good_z_sdss*good_z_petro]-z_sdss[good_z_sdss*good_z_petro]
    offset_dm_z_petro = U.mean_robust(dm_z_petro)
    dm_z_aper = z_auto[good_z_sdss*good_z_aper]-z_sdss[good_z_sdss*good_z_aper]
    offset_dm_z_aper = U.mean_robust(dm_z_aper) 
    # Saving info.
    line  = 'zSDSS  %.3f  %.3f '%(offset_dm_z_auto,offset_dm_z_petro)
    line +=' %.3f \n'%(offset_dm_z_aper)
    filename.write(line)    
    
    # Closing file.
    filename.close()
