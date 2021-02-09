__author__ = 'albertomolino'

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats/'
redshift = '05'
master_bpz_auto  = root2cats + 'master_STRIPE82_spz.z%s_cali_GOSMOS_auto.bpz'%(redshift)
master_bpz_petro = root2cats + 'master_STRIPE82_spz.z%s_cali_GOSMOS_petro.bpz'%(redshift)
master_bpz_aper  = root2cats + 'master_STRIPE82_spz.z%s_cali_GOSMOS_aper.bpz'%(redshift)
ng = len(U.get_data(master_bpz_auto,0))

######## ######### #########
#  New figures
######## ######### #########

# Precision as a function of z for several AB intervals.
# Version II

#base_z  = N.arange(0.,0.51,0.1)
base_z  = N.arange(0.0,0.25,0.05)
base_m  = N.arange(14.,21.,1.)
base_z2 = base_z[:-1]+((base_z[1]-base_z[0])/2.)
base_m2 = base_m[:-1]+((base_m[1]-base_m[0])/2.)
n_z = len(base_z)-1
n_m = len(base_m)-1
valor_auto = N.zeros((n_m,n_z),float)

alphas = [0.2,0.4,0.6,0.8,1.0]

#AUTO apertures
zb,zs,mo,ods = U.get_data(master_bpz_auto,(1,9,10,5))
dz = (zb-zs)/(1.+zs)

for ii in range(n_m):
    print ' '
    for jj in range(n_z):
        good  = N.greater_equal(mo,base_m[ii])
        good *= N.less_equal(mo,base_m[ii+1])
        good *= N.greater_equal(zs,base_z[jj])
        good *= N.less_equal(zs,base_z[jj+1])
        good *= N.greater_equal(ods,0.)
        if len(dz[good])>100:
            valor_auto[ii,jj] = U.std_mad(dz[good])
        else:
            valor_auto[ii,jj] = -1.

# Plot
plt.figure(1,figsize=(9,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
for ii in range(n_m):
    for jj in range(n_z):
        if valor_auto[ii,jj]>0:
           plt.scatter(base_z2[jj],base_m2[ii],s=4000,
                    c=valor_auto[ii,jj],marker=u's',
                    cmap=cm.PuOr,alpha=0.95,
                    vmin=0.01,vmax=0.04)

#cmap=cm.PuOr,alpha=0.85,
#cm.jet
cb = plt.colorbar(pad=0.,format='%.2f',ticks=[0.01,0.02,0.03,0.04,0.05,0.06])
cb.set_label('Photometric Redshift Precision',size=25,labelpad=10)
plt.grid()
plt.xlim(0.0,0.2)
plt.ylim(14,20)
#plt.legend(label4legend,loc='lower right',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('$r$',size=30,labelpad=4)
plt.xlabel('$z$',size=35,labelpad=-1)



# Precision as a function of z for several AB intervals.

base_m  = N.array([16,18,19,20])
base_z  = N.array([0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5])
base_z2 = base_z[:-1]+((base_z[1]-base_z[0])/2.)
base_m2 = base_m[:-1]+((base_m[1]-base_m[0])/2.)
n_z = len(base_z)
n_m = len(base_m)
valor_auto = N.zeros((n_m,n_z),float)

alphas = [0.2,0.4,0.6,0.8,1.0]

#AUTO apertures
zb,zs,mo = U.get_data(master_bpz_auto,(1,9,10))
dz = (zb-zs)/(1.+zs)
satur = N.greater(mo,14)

for ii in range(n_m):
    good_1 = N.less_equal(mo,base_m[ii])
    zs_r,dz_r = U.multicompress(good_1,(zs,dz))
    for jj in range(n_z):
        #good = N.greater_equal(zs_r,base_z[jj])
        good = N.less_equal(zs_r,base_z[jj])
        linea  = 'm<%.2f,'%(base_m[ii])
        linea += 'z<%.2f'%(base_z[jj])
        linea += ': %i'%(len(zs_r[good]))
        if len(zs_r[good])>100:
            valor_auto[ii,jj] = U.std_mad(dz_r[good])
        else:
            valor_auto[ii,jj] = -1.
        linea += ', dz/1+z: %.3f'%(valor_auto[ii,jj])
        print linea
        #pausa = raw_input('paused')

# Plot
plt.figure(1,figsize=(9,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
label4legend = ['R<16','R<18','R<19','R<20']
valor_auto[0,0]=0.012
for ii in range(n_m):
    sigma_z = valor_auto[ii,:]
    non_zero = N.greater(sigma_z,-1)
    base_redu,sigma_z = U.multicompress(non_zero,(base_z,sigma_z))
    plt.plot(base_redu,sigma_z,'-s',lw=10,alpha=alphas[ii],ms=12,color='#660033')
    #label4legend.append()
plt.grid()
plt.legend(label4legend,loc='lower right',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('$\sigma_{z}$',size=30,labelpad=-5)
plt.xlabel('$z$',size=35,labelpad=10)


# Precision as a function of AB for several z intervals
base_m  = N.array([15,16,17,18,19,20,21])
base_z  = N.array([0.1,0.3,0.4,0.5])
base_z2 = base_z[:-1]+((base_z[1]-base_z[0])/2.)
base_m2 = base_m[:-1]+((base_m[1]-base_m[0])/2.)
n_z = len(base_z)
n_m = len(base_m)
valor_auto = N.zeros((n_z,n_m),float)

alphas = [0.2,0.4,0.6,0.8,1.0]

#AUTO apertures
zb,zs,mo = U.get_data(master_bpz_auto,(1,9,10))
satu = N.greater(mo,14)
zb,zs,mo = U.multicompress(satu,(zb,zs,mo))
dz = (zb-zs)/(1.+zs)
for ii in range(n_z):
    good_1 = N.less_equal(zs,base_z[ii])
    zs_r,dz_r,mo_r = U.multicompress(good_1,(zs,dz,mo))
    for jj in range(n_m):
        #good = N.greater_equal(zs_r,base_z[jj])
        good = N.less_equal(mo_r,base_m[jj])
        linea  = 'm<%.2f,'%(base_z[ii])
        linea += 'z<%.2f'%(base_m[jj])
        linea += ': %i'%(len(mo_r[good]))
        if len(mo_r[good])>100:
            valor_auto[ii,jj] = U.std_mad(dz_r[good])
        else:
            valor_auto[ii,jj] = -1.
        linea += ', dz/1+z: %.3f'%(valor_auto[ii,jj])
        print linea
        #pausa = raw_input('paused')

# Plot
plt.figure(2,figsize=(9,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
label4legend = ['z<0.1','z<0.2','z<0.3','z<0.4','z<0.5']
#valor_auto[0,0]=0.012
for ii in range(n_z):
    sigma_z = valor_auto[ii,:]
    try:
       non_zero = N.greater(sigma_z,-1)
       base_redu,sigma_z = U.multicompress(non_zero,(base_z,sigma_z))
    except:
       base_redu = basem *1.
    plt.plot(base_redu,sigma_z,'-s',lw=10,alpha=alphas[ii],ms=12,color='#660033')
    #label4legend.append()
plt.grid()
plt.legend(label4legend,loc='lower right',fontsize=20)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('$\sigma_{z}$',size=30,labelpad=-5)
plt.xlabel('$AB$',size=35,labelpad=10)




"""
z<0.1
14-15,15-16,16-17,17-18,18-19,19-20
0.013,0.015,0.018,0.024,0.032,0.042
z<0.2
14-15,15-16,16-17,17-18,18-19,19-20
0.013,0.016,0.019,0.023,0.037,0.053
z<0.3
14-15,15-16,16-17,17-18,18-19,19-20
0.013,0.016,0.019,0.024,0.035,0.052
z<0.4
14-15,15-16,16-17,17-18,18-19,19-20
0.013,0.016,0.019,0.024,0.034,0.046
z<0.5
14-15,15-16,16-17,17-18,18-19,19-20
0.013,0.016,0.019,0.024,0.034,0.043
"""
###################################################
###################################################
###################################################

master_bpz_auto = root2cats + 'master_STRIPE82_spz.z%s_cali_GOSMOS_auto.bpz'%(redshift)
#master_bpz_auto  = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
#master_bpz_auto += 'splus_cats/BPZ_5_Filters/master_Stripe82_5bands_cali_GOSMOSeB11.bpz'

### General numbers
z_max=0.3
m_min=19.
m_max=20.
o_min=0.
t_min=0
t_max=100
a1=B.d_stats(master_bpz_auto,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min,tmin=t_min,tmax=t_max).nice().split('\n')[1]
#a2=B.d_stats(master_bpz_petro,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min,tmin=t_min,tmax=t_max).nice().split('\n')[1]
#a3=B.d_stats(master_bpz_aper,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min,tmin=t_min,tmax=t_max).nice().split('\n')[1]
print a1+' (%.2f)  AUTO'%(float(a1.split()[-1])/(1.*ng))
#print a2+' (%.2f)  PETRO'%(float(a2.split()[-1])/(1.*ng))
#print a3+' (%.2f)  APER'%(float(a3.split()[-1])/(1.*ng))


### Precision as a function of AB
m_min = 14.
m_max = 20.
delta_m = 0.5
base_m = N.arange(m_min,m_max+delta_m,delta_m)
base_m2 = base_m[:-1]+((base_m[1]-base_m[0])/2.)

#AUTO
zb,zs,mo = U.get_data(master_bpz_auto,(1,9,10))
dz = (zb-zs)/(1.+zs)
valor_auto = N.zeros(len(base_m)-1)
for ii in range(len(valor_auto)):
    good  = N.greater_equal(mo,base_m[ii])
    good *= N.less_equal(mo,base_m[ii+1])
    valor_auto[ii] = U.std_mad(dz[good])

#PETRO
zb,zs,mo = U.get_data(master_bpz_petro,(1,9,10))
dz = (zb-zs)/(1.+zs)
valor_petro = N.zeros(len(base_m)-1)
for ii in range(len(valor_petro)):
    good  = N.greater_equal(mo,base_m[ii])
    good *= N.less_equal(mo,base_m[ii+1])
    valor_petro[ii] = U.std_mad(dz[good])

#APER
zb,zs,mo = U.get_data(master_bpz_aper,(1,9,10))
dz = (zb-zs)/(1.+zs)
valor_aper = N.zeros(len(base_m)-1)
for ii in range(len(valor_aper)):
    good  = N.greater_equal(mo,base_m[ii])
    good *= N.less_equal(mo,base_m[ii+1])
    valor_aper[ii] = U.std_mad(dz[good])


plt.figure(1,figsize=(9,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
#plt.plot(base_m2,valor_aper,'-o',color='green',lw=8,alpha=0.5,ms=20)
#plt.plot(base_m2,valor_petro,'-o',color='red',lw=8,alpha=0.5,ms=20)
plt.plot(base_m2,valor_auto,'-o',color='red',lw=8,alpha=0.5,ms=20)
plt.grid()
plt.legend(['aper','petro','auto'],loc='best',fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('$\sigma_{z}$',size=30,labelpad=-5)
plt.xlabel('$r$ [AB]',size=25,labelpad=10)


### Precision as a function of z
z_min = 0.005
z_max = 0.2
delta_z = 0.005
base_z = N.arange(z_min,z_max+delta_z,delta_z)
base_z2 = base_z[:-1]+((base_z[1]-base_z[0])/2.)

#AUTO
zb,zs,mo = U.get_data(master_bpz_auto,(1,9,10))
dz = (zb-zs)/(1.+zs)
valor_auto = N.zeros(len(base_z)-1)
for ii in range(len(valor_auto)):
    good  = N.greater_equal(zs,base_z[ii])
    good *= N.less_equal(zs,base_z[ii+1])
    valor_auto[ii] = U.std_mad(dz[good])

#PETRO
zb,zs,mo = U.get_data(master_bpz_petro,(1,9,10))
dz = (zb-zs)/(1.+zs)
valor_petro = N.zeros(len(base_z)-1)
for ii in range(len(valor_petro)):
    good  = N.greater_equal(zs,base_z[ii])
    good *= N.less_equal(zs,base_z[ii+1])
    valor_petro[ii] = U.std_mad(dz[good])

#APER
zb,zs,mo = U.get_data(master_bpz_aper,(1,9,10))
dz = (zb-zs)/(1.+zs)
valor_aper = N.zeros(len(base_z)-1)
for ii in range(len(valor_aper)):
    good  = N.greater_equal(zs,base_z[ii])
    good *= N.less_equal(zs,base_z[ii+1])
    valor_aper[ii] = U.std_mad(dz[good])


plt.figure(1,figsize=(9,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
plt.plot(base_z2,valor_aper,'-o',color='green',lw=8,alpha=0.5,ms=20)
plt.plot(base_z2,valor_petro,'-o',color='red',lw=8,alpha=0.5,ms=20)
plt.plot(base_z2,valor_auto,'-o',color='blue',lw=8,alpha=0.5,ms=20)
plt.grid()
plt.legend(['aper','petro','auto'],loc='best',fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('$\sigma_{z}$',size=30,labelpad=-5)
plt.xlabel('$z_{s}$',size=25,labelpad=10)


### Precision as a function of Odds
base_o = N.array([0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.55,
                  0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,0.99])
base_o2 = base_o[:-1]+((base_o[1]-base_o[0])/2.)

#AUTO
zb,zs,odds = U.get_data(master_bpz_auto,(1,9,5))
dz = (zb-zs)/(1.+zs)
#g  = N.less_equal(zs,0.15)
#dz = dz[g]
valor_auto = N.zeros(len(base_o)-1)
for ii in range(len(valor_auto)):
    good  = N.greater_equal(odds,base_o[ii])
    good *= N.less_equal(odds,base_o[ii+1])
    valor_auto[ii] = U.std_mad(dz[good])

#PETRO
zb,zs,odds = U.get_data(master_bpz_petro,(1,9,5))
dz = (zb-zs)/(1.+zs)
valor_petro = N.zeros(len(base_o)-1)
for ii in range(len(valor_petro)):
    good  = N.greater_equal(odds,base_o[ii])
    good *= N.less_equal(odds,base_o[ii+1])
    valor_petro[ii] = U.std_mad(dz[good])

#APER
zb,zs,odds = U.get_data(master_bpz_aper,(1,9,5))
dz = (zb-zs)/(1.+zs)
valor_aper = N.zeros(len(base_o)-1)
for ii in range(len(valor_aper)):
    good  = N.greater_equal(odds,base_o[ii])
    good *= N.less_equal(odds,base_o[ii+1])
    valor_aper[ii] = U.std_mad(dz[good])


plt.figure(1,figsize=(9,8),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
plt.plot(base_o2,valor_aper,'-o',color='green',lw=8,alpha=0.5,ms=20)
plt.plot(base_o2,valor_petro,'-o',color='red',lw=8,alpha=0.5,ms=20)
plt.plot(base_o2,valor_auto,'-o',color='blue',lw=8,alpha=0.5,ms=20)
plt.grid()
plt.legend(['aper','petro','auto'],loc='best',fontsize=30)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.ylabel('$\sigma_{z}$',size=30,labelpad=-5)
plt.xlabel('$Odds$',size=25,labelpad=10)

### Completeness in z as a function of Odds. pepe


### Comparison among PDZ and dz/1+z





"""
root ='/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/splus_cats_NGSL/'
b0 = root + 'COSMOSeB11new_recal/PriorSM/ALL_BPZ_OLD/master.bpz'
b1 = root + 'May18_recal/NoSEDrecal/masterBPZnoSEDrecal.bpz'
b2 = root + 'May18_recal/masterALL.bpz'

z_min=0.;z_max=1.;m_min=14.;m_max=22.;o_min=0.0;t_min=0;t_max=100.5;chimax=30;ccut=0.05;a0=B.d_stats(b0,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min,tmin=t_min,tmax=t_max,cut=ccut,chi2max=chimax).nice().split('\n')[1];a1=B.d_stats(b1,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min,tmin=t_min,tmax=t_max,cut=ccut,chi2max=chimax).nice().split('\n')[1];a2=B.d_stats(b2,zmin=z_min,zmax=z_max,mmin=m_min,mmax=m_max,omin=o_min,tmin=t_min,tmax=t_max,chi2max=chimax,cut=ccut).nice().split('\n')[1];print a0+'ORICAL';print a1+'NEWCAL NOSEDrecal';print a2+'NEWCAL SEDrecal'

"""