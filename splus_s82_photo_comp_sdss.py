__author__ = 'albertomolino'

# /Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/
import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root += 'data_quality/photometry/SDSSSPLUS/'
u_cat = root + 'master_SPLUS_STRIPE82_to_SDSS_Ivezick_stars_Uband.cat'
griz_cat = root + 'master_SPLUS_STRIPE82_to_SDSS_Ivezick_stars_griz.cat'


s2nu,u_sp,u_sd = U.get_data(u_cat,(2,3,4))
delta_u = (u_sp-u_sd)
good_u = N.greater_equal(s2nu,5.) * N.less_equal(abs(delta_u),2.)
good_u *= N.less_equal(u_sd,22.)
u_sp,u_sd,delta_u = U.multicompress(good_u,(u_sp,u_sd,delta_u))

s2n,g_sp,r_sp,i_sp,z_sp = U.get_data(griz_cat,(2,3,4,5,6)) # SPLUS
g_sd,r_sd,i_sd,z_sd = U.get_data(griz_cat,(7,8,9,10))  # SDSS
delta_g = (g_sp-g_sd)
delta_r = (r_sp-r_sd)
delta_i = (i_sp-i_sd)
delta_z = (z_sp-z_sd)

good_griz = N.greater_equal(s2n,5.) * N.less_equal(abs(delta_g),2.)
good_griz *= N.less_equal(abs(delta_r),2.) * N.less_equal(abs(delta_i),2.)
good_griz *= N.less_equal(abs(delta_z),2.)
g_sp,r_sp,i_sp,z_sp = U.multicompress(good_griz,(g_sp,r_sp,i_sp,z_sp))
g_sd,r_sd,i_sd,z_sd = U.multicompress(good_griz,(g_sd,r_sd,i_sd,z_sd))
delta_g,delta_r = U.multicompress(good_griz,(delta_g,delta_r))
delta_i,delta_z = U.multicompress(good_griz,(delta_i,delta_z))

base_color = N.arange(-2,2.,0.02)
res = 5
#splus_colors = list(cm.jet(N.linspace(0, 1, 12)))
#splus_colors = list(cm.rainbow(N.linspace(0, 1, 5)))
y_min = -1.49
y_max =  1.49
y_min = -0.99
y_max =  0.99
x_min =  14
x_max =  21

## U
#plt.figure(1,figsize = (14,6),dpi=70, facecolor='w', edgecolor='k')
plt.figure(1,figsize = (14,4),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
uno = plt.axes([.075,.12,.725,.87])
plt.plot(u_sd[::res],delta_u[::res]-(U.mean_robust(delta_u)),'o',
         ms=5,alpha=0.1,color='grey')
plt.xlim(x_min,x_max)
plt.ylim(y_min,y_max)
plt.grid()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel('$m^{sp}_{u}$ - $m^{sd}_{u}$',size=25,labelpad=-1)
#plt.xlabel('$m^{sd}_{u}$',size=25,labelpad=1)
plt.xlabel('$m^{u,sp}$',size=23,labelpad=-10)

dos = plt.axes([.8,.12,.15,.87])
a1,a2,a3 = plt.hist(delta_u[::res]-(U.mean_robust(delta_u)),
                    base_color,orientation='horizontal',color='grey',
                    alpha=0.5,normed=True)
plt.grid()
plt.ylim(y_min,y_max)
plt.xticks([])
plt.yticks([])
plt.savefig('/Users/albertomolino/Desktop/dm_u.png',dpi=100)

## G
plt.figure(2,figsize = (14,4),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
uno = plt.axes([.075,.12,.725,.87])
plt.plot(g_sd[::res],delta_g[::res]-(U.mean_robust(delta_g[g_sp<19])),'o',
         ms=5,alpha=0.1,color='grey')
plt.xlim(x_min,x_max)
plt.ylim(y_min,y_max)
plt.grid()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel('$m^{sp}_{g}$ - $m^{sd}_{g}$',size=25,labelpad=-1)
plt.xlabel('$m^{g,sd}$',size=23,labelpad=-10)

dos = plt.axes([.8,.12,.15,.87])
a1,a2,a3 = plt.hist(delta_g[::res]-(U.mean_robust(delta_g[g_sp<19])),
                    base_color,orientation='horizontal',color='grey'
                    ,alpha=0.5,normed=True)
plt.grid()
plt.ylim(y_min,y_max)
plt.xticks([])
plt.yticks([])
plt.savefig('/Users/albertomolino/Desktop/dm_g.png',dpi=100)

## R
plt.figure(3,figsize = (14,4),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
uno = plt.axes([.075,.12,.725,.87])
plt.plot(r_sd[::res],delta_r[::res]-(U.mean_robust(delta_r[r_sp<19])),'o',
         ms=5,alpha=0.1,color='grey')
plt.xlim(x_min,x_max)
plt.ylim(y_min,y_max)
plt.grid()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel('$m^{sp}_{r}$ - $m^{sd}_{r}$',size=25,labelpad=-1)
plt.xlabel('$m^{r,sd}$',size=23,labelpad=-10)

dos = plt.axes([.8,.12,.15,.87])
a1,a2,a3 = plt.hist(delta_r[::res]-(U.mean_robust(delta_r[r_sp<19])),
                    base_color,orientation='horizontal',color='grey',
                    alpha=0.5,normed=True)
plt.grid()
plt.ylim(y_min,y_max)
plt.xticks([])
plt.yticks([])
plt.savefig('/Users/albertomolino/Desktop/dm_r.png',dpi=100)

## I
plt.figure(4,figsize = (14,4),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
uno = plt.axes([.075,.12,.725,.87])
plt.plot(i_sd[::res],delta_i[::res]-(U.mean_robust(delta_i[i_sp<19])),'o',
         ms=5,alpha=0.1,color='grey')
plt.xlim(x_min,x_max)
plt.ylim(y_min,y_max)
plt.grid()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel('$m^{sp}_{i}$ - $m^{sd}_{i}$',size=25,labelpad=-1)
plt.xlabel('$m^{i,sd}$',size=23,labelpad=-10)

dos = plt.axes([.8,.12,.15,.87])
a1,a2,a3 = plt.hist(delta_i[::res]-(U.mean_robust(delta_i[i_sp<19])),
                    base_color,orientation='horizontal',color='grey',
                    alpha=0.5,normed=True)
plt.grid()
plt.ylim(y_min,y_max)
plt.xticks([])
plt.yticks([])
plt.savefig('/Users/albertomolino/Desktop/dm_i.png',dpi=100)

## Z
plt.figure(5,figsize = (14,4),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
uno = plt.axes([.075,.12,.725,.87])
plt.plot(z_sd[::res],delta_z[::res]-(U.mean_robust(delta_z[z_sp<19])),'o',
         ms=5,alpha=0.1,color='grey')
plt.xlim(x_min,x_max)
plt.ylim(y_min,y_max)
plt.grid()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.ylabel('$m^{sp}_{z}$ - $m^{sd}_{z}$',size=25,labelpad=-1)
plt.xlabel('$m^{z,sd}$',size=23,labelpad=-10)

dos = plt.axes([.8,.12,.15,.87])
a1,a2,a3 = plt.hist(delta_z[::res]-(U.mean_robust(delta_z[z_sp<19])),
                    base_color,orientation='horizontal',color='grey',
                    alpha=0.5,normed=True)
plt.grid()
plt.ylim(y_min,y_max)
plt.xticks([])
plt.yticks([])
plt.savefig('/Users/albertomolino/Desktop/dm_z.png',dpi=100)



print 'std_mad_U',U.std_mad(delta_u[u_sp<19])
print 'std_mad_G',U.std_mad(delta_g[g_sp<19])
print 'std_mad_R',U.std_mad(delta_r[r_sp<19])
print 'std_mad_I',U.std_mad(delta_i[i_sp<19])
print 'std_mad_Z',U.std_mad(delta_z[z_sp<19])

"""
std_mad_U 0.059304
std_mad_G 0.0474432
std_mad_R 0.0266868
std_mad_I 0.0474432
std_mad_Z 0.0311346
"""
