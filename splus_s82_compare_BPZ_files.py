__author__ = 'albertomolino'

"""
This routine shows the performance of several BPZ files,
derived using different SEDs and Priors.

"""

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import bpz_tools as B

#Roots
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/'

# Files
#bpz_0 = root2cats + 'bpz/master.STRIPE82_Photometry.m21.bpz'
#bpz_1 = root2cats + 'newSED_v1/SEDrecal_sorted/master.STRIPE82_Photometry.m21.bpz'
#bpz_2 = root2cats + 'newSED_v1/new_SEDsorted_recal_all/master.STRIPE82_Photometry.m21.recal_all.bpz'
#bpz_3 = root2cats + 'newSED_v1/new_SEDsorted_recal_2reds/master.STRIPE82_Photometry.m21.recal_2reds.bpz'
#bpz_4 = root2cats + 'newSED_v1/new_SEDsorted_no_further_recal/master.STRIPE82_Photometry.m21.nofurhtercal.bpz'
#bpz_5 = root2cats + 'COSMOSeB11new/master.STRIPE82_Photometry.m21_COSMOSeB11new.bpz'
bpz_6 = root2cats + 'COSMOSeB11new_recal/master.STRIPE82_Photometry.m21_COSMOSeB11new_recal.bpz'
bpz_7 = root2cats + 'COSMOSeB11new_recal/master.STRIPE82_Photometry.m21_COSMOSeB11new_recal_redu.bpz'
#Selection input value
z_min=0.0
z_max=0.5
chi2_max=30
m_min=14.
m_max=21.
o_min=0.05
t_min=0.
t_max=50.

#Executing...
#a0 = B.d_stats(bpz_0,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).nice().split('\n')[1]
#a1 = B.d_stats(bpz_1,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).nice().split('\n')[1]
#a2 = B.d_stats(bpz_2,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).nice().split('\n')[1]
#a3 = B.d_stats(bpz_3,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).nice().split('\n')[1]
#a4 = B.d_stats(bpz_4,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).nice().split('\n')[1]
#a5 = B.d_stats(bpz_5,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).nice().split('\n')[1]
a6 = B.d_stats(bpz_6,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).nice().split('\n')[1]
a7 = B.d_stats(bpz_7,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).nice().split('\n')[1]

#print a0+' v0'
#print a1+' v1'
#print a2+' v2all'
#print a3+' v2reds'
#print a4+' v2nocals'
#print a5+' COSMOSeB11new'
print a6+' COSMOSeB11new_recal'
print a7+' COSMOSeB11new_recal_redu'

#a0 = B.d_stats(bpz_0,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).types().split('\n')
#a1 = B.d_stats(bpz_1,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).types().split('\n')
#a5 = B.d_stats(bpz_5,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).types().split('\n')
a6 = B.d_stats(bpz_6,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).types().split('\n')
a7 = B.d_stats(bpz_7,chi2max=chi2_max,tmin=t_min,tmax=t_max,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min).types().split('\n')


for ii in range(16):
    if ii == 0: print ''#'  v0                   v1                         v2all              v2reds'
    #else: print a0[ii],'  ',a1[ii],'  ',a2[ii],'  ',a3[ii],'  ',a4[ii],'  ',a5[ii],'  ',a6[ii]
    else: print a6[ii],'  ',a7[ii],#'  ',a5[ii],'  ',a6[ii]

