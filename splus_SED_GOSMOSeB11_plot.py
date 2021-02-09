__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt
from matplotlib import cm

# Paths
root_to_bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root_to_seds = root_to_bpz+'SED/'
sed = U.get_str(root_to_seds+'COSMOSeB11new_recal.list',0)
n_seds = len(sed)-2

def lookcloser(vector, value):
    dim = len(vector)
    try:
       if vector[0] >= value:    pos = 0
       elif vector[1] >= value:  pos = 1
       elif vector[-1:] < value: pos = dim-1
       else:
         for ii in range(len(vector)):
           if vector[ii-2] < value < vector[ii]:
               pos = ii-1
    except:
          print 'ID not found!'
          pos  = -1
    return pos

seds_colors = list(cm.jet(N.linspace(0, 1, n_seds)))
seds_colors = seds_colors[::-1]

# Plots
plt.figure(10, figsize = (11,9),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
for ii in range(n_seds):
    x,y = U.get_data(root_to_seds+sed[ii],(0,1))
    pepe = lookcloser(x,4000.)
    plt.loglog(x,y/y[pepe-1],'-',linewidth=4.5,alpha=0.5,color=seds_colors[ii])
    plt.xlim(800,20000)
    plt.ylim(0.001,10)

for ii in range(n_seds):
    x,y = U.get_data(root_to_seds+sed[ii],(0,1))
    pepe = lookcloser(x,4000.)
    plt.loglog(x,y/y[pepe-1],'-',linewidth=1.5,alpha=0.99,color=seds_colors[ii])
    plt.xlim(800,20000)
    plt.ylim(0.001,10)

plt.xticks(fontsize=19)
plt.yticks(fontsize=19)
plt.xlabel('Wavelength [$\AA$]',size=21,labelpad=-3)
plt.ylabel('Flux  [erg/s/cm$^{2}$/$\AA$]',size=21,labelpad=-3)

leg = ['ELL','ELL','ELL','ELL','ES0','Sa','Sa','Sa','Sb','Sb','Sbc','Scd','SB','SB']
plt.legend(leg,loc='best',fontsize=17,frameon=False)


"""
res=2
plt.figure(2)
plt.clf()
plt.scatter(mo[::res],t[::res],s=100,
c=chi[::res],marker=u's',cmap=cm.PuOr,
alpha=0.1,vmin=0.,vmax=2.)
cb = plt.colorbar(pad=0.,format='%.2f')
plt.xlim(14,21)
plt.ylim(0.,20)
plt.grid()

-
colores[:,0] =(0.50,0.00,0.00)
colores[:,1] =(0.65,0.00,0.00)
colores[:,2] =(0.75,0.00,0.00)
colores[:,3] =(0.85,0.00,0.00)
colores[:,4] =(0.90,0.00,0.00)
colores[:,5] =(1.00,0.00,0.00)
colores[:,6] =(0.80,0.05,0.00)
colores[:,7] =(0.80,0.25,0.00) ##
colores[:,8] =(0.85,0.25,0.00)
colores[:,9] =(0.75,0.50,0.00)
colores[:,10]=(0.85,0.50,0.00) ##
colores[:,11]=(0.85,0.65,0.00)
colores[:,12]=(0.75,0.75,0.00)
colores[:,13]=(0.00,0.50,0.00)
colores[:,14]=(0.00,0.50,0.00)
colores[:,15]=(0.00,0.65,1.00)
colores[:,16]=(0.00,0.50,1.00)
colores[:,17]=(0.00,0.25,1.00)
colores[:,18]=(0.00,0.20,1.00)
colores[:,19]=(0.00,0.00,1.00)

plt.legend(leg,loc='upper right', bbox_to_anchor=(0.5, 0.5),fontsize=15,frameon=False,labelspacing=0.2,framealpha=0.3)

#Color definition.
colores = N.zeros((3,15),float)
#colores[:,0] =(0.50,0.00,0.00)
#colores[:,1] =(0.65,0.00,0.00)
#colores[:,2] =(0.75,0.00,0.00)
#colores[:,3] =(0.85,0.00,0.00)
#colores[:,4] =(0.90,0.00,0.00)
colores[:,0] =(1.00,0.00,0.00)
colores[:,1] =(0.80,0.05,0.00)
colores[:,2] =(0.80,0.25,0.00) ##
colores[:,3] =(0.85,0.25,0.00)
colores[:,4] =(0.75,0.50,0.00)
colores[:,5]=(0.85,0.50,0.00) ##
colores[:,6]=(0.85,0.65,0.00)
colores[:,7]=(0.75,0.75,0.00)
colores[:,8]=(0.00,0.50,0.00)
colores[:,9]=(0.00,0.50,0.00)
colores[:,10]=(0.00,0.65,1.00)
colores[:,11]=(0.00,0.50,1.00)
colores[:,12]=(0.00,0.25,1.00)
colores[:,13]=(0.00,0.20,1.00)
colores[:,14]=(0.00,0.00,1.00)

"""