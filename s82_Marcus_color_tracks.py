import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import bpz_tools as B
import numpy as N

# This is the root with the original files.
root = '/Users/albertomolino/doctorado/articulos/SPLUS/StarGalaxy/AB/'

"""
Ell3_A_0.SPLUS_gSDSS.AB.txt
Ell3_A_0.SPLUS_iSDSS.AB.txt	
Ell3_A_0.SPLUS_rSDSS.AB.txt	
Ell3_A_0.SPLUS_zSDSS.AB.txt	
Sa_A_1.SPLUS_gSDSS.AB.txt
Sa_A_1.SPLUS_iSDSS.AB.txt
Sa_A_1.SPLUS_rSDSS.AB.txt
Sa_A_1.SPLUS_zSDSS.AB.txt
"""
zmax = 2.

#sed = 'Ell3_A_0'
sed = 'Sa_A_1'

filtro1 = 'rSDSS'
filtro2 = 'iSDSS'

# filtro1 = 'gSDSS'
# filtro2 = 'rSDSS'
# filtro1 = 'iSDSS'
# filtro2 = 'zSDSS'

sed_ab_filter_1 = root + '%s.SPLUS_%s.AB.txt'%(sed,filtro1)
sed_ab_filter_2 = root + '%s.SPLUS_%s.AB.txt'%(sed,filtro2)

z_AB_1,f_AB_1 = U.get_data(sed_ab_filter_1,(0,1)) 
z_AB_2,f_AB_2 = U.get_data(sed_ab_filter_2,(0,1))

good_z_range = N.less_equal(z_AB_1,zmax)
good_z_range *= N.greater(f_AB_1,0) * N.greater(f_AB_2,0)
f_AB_redu_1 = f_AB_1[good_z_range]
z_AB_redu_1 = z_AB_1[good_z_range]
f_AB_redu_2 = f_AB_2[good_z_range]
z_AB_redu_2 = z_AB_2[good_z_range]

mag_AB_redu_1 = B.flux2mag(abs(f_AB_redu_1))
mag_AB_redu_2 = B.flux2mag(abs(f_AB_redu_2))

outfilename = root + 'tracks_sed%s_0.0z%.1f_%s%s.txt'%(sed,zmax,filtro1,filtro2)
U.put_data(outfilename,(z_AB_redu_2,mag_AB_redu_1-mag_AB_redu_2),'z %s-%s'%(filtro1,filtro2))


"""
z,gr_ell = U.get_data('/Users/albertomolino/doctorado/articulos/SPLUS/StarGalaxy/AB/tracks_sedEll3_A_0_0.0z0.5_gSDSSrSDSS.txt',(0,1));z,rz_ell = U.get_data('/Users/albertomolino/doctorado/articulos/SPLUS/StarGalaxy/AB/tracks_sedEll3_A_0_0.0z0.5_rSDSSzSDSS.txt',(0,1));z,gr_Sa = U.get_data('/Users/albertomolino/doctorado/articulos/SPLUS/StarGalaxy/AB/tracks_sedSa_A_1_0.0z0.5_gSDSSrSDSS.txt',(0,1)); z,rz_Sa = U.get_data('/Users/albertomolino/doctorado/articulos/SPLUS/StarGalaxy/AB/tracks_sedEll3_A_0_0.0z0.5_rSDSSzSDSS.txt',(0,1)); plt.plot(gr_Sa[::10],rz_Sa[::10],'-bo');plt.plot(gr_ell[::10],rz_ell[::10],'-ro')

z,gr_ell = U.get_data('/Users/albertomolino/doctorado/articulos/SPLUS/StarGalaxy/AB/tracks_sedEll3_A_0_0.0z2.0_gSDSSrSDSS.txt',(0,1));z,rz_ell = U.get_data('/Users/albertomolino/doctorado/articulos/SPLUS/StarGalaxy/AB/tracks_sedEll3_A_0_0.0z2.0_rSDSSzSDSS.txt',(0,1));z,gr_Sa = U.get_data('/Users/albertomolino/doctorado/articulos/SPLUS/StarGalaxy/AB/tracks_sedSa_A_1_0.0z2.0_gSDSSrSDSS.txt',(0,1)); z,rz_Sa = U.get_data('/Users/albertomolino/doctorado/articulos/SPLUS/StarGalaxy/AB/tracks_sedSa_A_1_0.0z2.0_rSDSSzSDSS.txt',(0,1)); plt.plot(gr_Sa[::10],rz_Sa[::10],'-bo');plt.plot(gr_ell[::10],rz_ell[::10],'-ro')


plt.xlabel('$g-r$',size=25,labelpad=1)
plt.ylabel('$r-z$',size=25,labelpad=1)
plt.grid()
plt.legend(['early','late'],loc='best',fontsize=30,numpoints=1)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.title('0.0$<z<$2.0',size=20)
plt.savefig('/Users/albertomolino/doctorado/articulos/SPLUS/StarGalaxy/AB/tracks.png',dpi=90)

"""
