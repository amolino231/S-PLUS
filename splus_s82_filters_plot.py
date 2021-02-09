__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')

import useful as U
import numpy as N
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

#ruta = '/Users/albertomolino/codigos/bpz-1.99.2/FILTER/SPLUS_July2017/'
ruta = '/Users/albertomolino/codigos/bpz-1.99.2/FILTER/SPLUS_September2018/'
filtros = U.get_str(ruta+'SPLUS_201809.list',0)

colores = N.zeros((3,12),float)
colores[:,0]=(0.00,0.00,1.00)
colores[:,1]=(0.00,0.25,1.00)
colores[:,2]=(0.00,0.65,1.00)
colores[:,3]=(0.00,0.50,0.00)
colores[:,4]=(0.85,0.65,0.00)
colores[:,5]=(0.75,0.50,0.00)
colores[:,6]=(0.80,0.25,0.00)
colores[:,7]=(1.00,0.00,0.00)
colores[:,8]=(0.85,0.00,0.00)
colores[:,9]=(0.65,0.00,0.00)
colores[:,10]=(0.35,0.00,0.00)
colores[:,11]=(0.25,0.00,0.00)
NB=N.array([1,2,3,4,6,8,10])

plt.figure(1,figsize = (11.5,8.7),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
#plt.subplot(211)
for ii in range(12):
     x,y = U.get_data(ruta+filtros[ii],(0,1))
     plt.plot(x*1000.,y,'-',lw=10,color=colores[:,ii],alpha=0.5)


for ii in range(12):
     x,y = U.get_data(ruta+filtros[ii],(0,1))
     #plt.plot(x*1.,y,'--',lw=2,color=colores[:,ii],alpha=0.7)
     if ii in NB: plt.fill_between(x*1.,y,0,alpha=0.4,color=colores[:,ii],lw=4)
     else: plt.fill_between(x*1.,y,0,alpha=0.2,color=colores[:,ii],lw=2)

plt.xlabel('Wavelength [$\AA$]',size=30,labelpad=5)
plt.ylabel('Throughput',size=30,labelpad=10)
plt.grid()

lab = ['uJAVA','J0378','J0395','J0410','J0430','gSDSS',
          'J0515','rSDSS','J0660','iSDSS','J0861','zSDSS']

plt.legend(lab,loc='upper right',fontsize=21,frameon=False)
# plt.xlim(3000,14500)
plt.xlim(3000,10500)
plt.ylim(0.001,0.999)
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)


#### Separated in two windows
ruta = '/Users/albertomolino/codigos/bpz-1.99.2/FILTER/SPLUS_July2017/'
filtros = U.get_str(ruta+'SPLUS.list',0)

bbs = [0,5,7,9,11]
nbs = [1,2,3,4,6,8,10]

splus_colors = list(cm.jet(N.linspace(0, 1, 12)))

colores_BB = N.zeros((3,5),float)
colores_NB = N.zeros((3,7),float)
colores_BB[:,0]=(0.00,0.00,1.00)
colores_BB[:,1]=(0.00,0.50,0.00)
colores_BB[:,2]=(0.85,0.65,0.00)
colores_BB[:,3]=(1.00,0.00,0.00)
colores_BB[:,4]=(0.35,0.00,0.00)

colores_NB[:,0]=(0.00,0.25,1.00)
colores_NB[:,1]=(0.00,0.65,1.00)
colores_NB[:,2]=(0.00,0.50,0.00)
colores_NB[:,3]=(0.85,0.65,0.00)
colores_NB[:,4]=(0.80,0.25,0.00)
colores_NB[:,5]=(0.85,0.00,0.00)
colores_NB[:,6]=(0.25,0.00,0.00)


plt.figure(11,figsize = (11.5,8.7),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
plt.subplot(211)
#broad-bands
for ii in range(5):
     x,y = U.get_data(ruta+filtros[bbs[ii]],(0,1))
     color_index = bbs[ii]
     plt.plot(x*1000.,y,'-',lw=10,color=splus_colors[color_index],alpha=0.9)

for ii in range(5):
     x,y = U.get_data(ruta+filtros[bbs[ii]],(0,1))
     color_index = bbs[ii]
     plt.plot(x*1.,y,'-',lw=2,color=splus_colors[color_index],alpha=0.9)
     plt.fill_between(x*1.,y,0,alpha=0.4,color=splus_colors[color_index],lw=4)

plt.ylabel('Throughput',size=28,labelpad=10)
#plt.grid()

lab = ['uJAVA','gSDSS','rSDSS','iSDSS','zSDSS']
plt.legend(lab,loc='upper right',fontsize=23,frameon=False)

plt.xlim(3000,11900)
plt.ylim(0.01,0.79)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

plt.subplot(212)
for ii in range(7):
     x,y = U.get_data(ruta+filtros[nbs[ii]],(0,1))
     color_index = nbs[ii]
     plt.plot(x*1000.,y,'-',lw=10,color=splus_colors[color_index],alpha=0.9)

for ii in range(7):
     x,y = U.get_data(ruta+filtros[nbs[ii]],(0,1))
     color_index = nbs[ii]
     plt.plot(x*1.,y,'-',lw=3,color=splus_colors[color_index],alpha=0.9)
     plt.fill_between(x*1.,y,0,alpha=0.4,color=splus_colors[color_index],lw=4)


lab = ['J0378','J0395','J0410','J0430',
          'J0515','J0660','J0861']

plt.legend(lab,loc='upper right',fontsize=23,frameon=False)
# plt.legend(lab2,loc='upper right',fontsize=16)
plt.xlabel('Wavelength [$\AA$]',size=28,labelpad=5)
plt.ylabel('Throughput',size=28,labelpad=10)
#plt.grid()
plt.xlim(3000,11900)
plt.ylim(0.01,0.79)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


###### single frame

bbs = [0,5,7,9,11]
nbs = [1,2,3,4,6,8,10]

plt.figure(12,figsize = (11.5,8.7),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
#broad-bands
for ii in range(12):
     x,y = U.get_data(ruta+filtros[ii],(0,1))
     #plt.plot(x*1000.,y,'-',lw=10,color=splus_colors[ii],alpha=0.9)
     plt.plot(x*1000.,y,'-',lw=10,color=colores[:,ii],alpha=0.9)

for ii in range(12):
     x,y = U.get_data(ruta+filtros[ii],(0,1))
     if ii in bbs:
        plt.plot(x*1.,y,'-',lw=2,color=colores[:,ii],alpha=0.8)
        plt.fill_between(x*1.,y,0,alpha=0.3,color=colores[:,ii],lw=4)
        #plt.plot(x*1.,y,'-',lw=2,color=splus_colors[ii],alpha=0.8)
        #plt.fill_between(x*1.,y,0,alpha=0.3,color=splus_colors[ii],lw=4)
     else:
        #plt.plot(x*1.,y,'-',lw=2,color=splus_colors[ii],alpha=0.9)
        #plt.plot(x*1.,y,'-',lw=2,color='black',alpha=0.5)
        #plt.fill_between(x*1.,y,0,alpha=0.8,color=splus_colors[ii],lw=4)
        plt.plot(x*1.,y,'-',lw=2,color='black',alpha=0.5)
        plt.fill_between(x*1.,y,0,alpha=0.8,color=colores[:,ii],lw=4)


#lab = ['uJAVA','J0378','J0395','J0410','J0430','gSDSS',
#          'J0515','rSDSS','J0660','iSDSS','J0861','zSDSS']
lab = ['uJAVA','J0378','J0395','J0410','J0430','    g',
          'J0515','    r','J0660','    i','J0861','    z']
plt.legend(lab,loc='upper right',fontsize=23,frameon=False)

plt.xlim(3000,11900)
plt.ylim(0.01,0.99)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.xlabel('Wavelength [$\AA$]',size=28,labelpad=5)
plt.ylabel('Throughput',size=28,labelpad=10)
















"""
import useful as U
import numpy as N
import matplotlib.pyplot as plt

# ruta = '/Users/albertomolino/codigos/bpz-1.99.2/FILTER/'
ruta = '/Users/albertomolino/codigos/bpz-1.99.2/FILTER/SPLUS_groups/'
filtros = U.get_str(ruta+'splusBBs.list',0)

colores = N.zeros((3,5),float)
colores[:,0]=(0.00,0.00,1.00)
colores[:,1]=(0.00,0.50,0.00)
colores[:,2]=(0.75,0.50,0.00)
colores[:,3]=(1.00,0.00,0.00)
colores[:,4]=(0.25,0.00,0.00)

colores2 = N.zeros((3,19),float)
colores2[:,0] =(0.00,0.00,1.00) ## 378
colores2[:,1] =(0.00,0.15,1.00) ## 395
colores2[:,2] =(0.00,0.25,1.00) ## 410
colores2[:,3] =(0.00,0.35,1.00) ## 430
colores2[:,4] =(0.00,0.45,1.00) #  446
colores2[:,5] =(0.00,0.55,0.00) #  463
colores2[:,6] =(0.00,0.65,0.00) #  480
colores2[:,7] =(0.00,0.75,0.00) #  497
colores2[:,8] =(0.00,0.50,0.00) ## 515
colores2[:,9] =(0.00,0.70,0.00) #  533
colores2[:,10]=(0.50,0.50,0.00) #  553
colores2[:,11]=(0.75,0.75,0.00) #  573
colores2[:,12]=(0.79,0.75,0.00) #  593
colores2[:,13]=(0.85,0.65,0.00) #  613
colores2[:,14]=(0.85,0.50,0.00) #  630
colores2[:,15]=(0.80,0.50,0.00) ## 660
colores2[:,16]=(0.75,0.50,0.00) #  676
colores2[:,17]=(0.85,0.25,0.00) #  691
colores2[:,18]=(0.85,0.00,0.00) ## 861
NB=N.array([0,1,2,3,8,15,18])


plt.figure(1,figsize = (11.5,8.7),dpi=70, facecolor='w', edgecolor='k')
plt.clf()
plt.subplot(211)
for ii in range(5):
     x,y = U.get_data(filtros[ii],(0,1))
     plt.plot(x*1000.,y,'-',lw=10,color=colores[:,ii],alpha=0.5)

for ii in range(5):
     x,y = U.get_data(filtros[ii],(0,1))
     plt.plot(x*10.,y,'-',lw=2,color=colores[:,ii],alpha=0.5)
     plt.fill_between(x*10.,y,0,alpha=0.3,color=colores[:,ii],lw=4)

plt.ylabel('Throughput',size=28,labelpad=10)
plt.grid()

lab = ['uJAVA','gSDSS','rSDSS','iSDSS','zSDSS']
# plt.legend(lab,loc='upper right',fontsize=32)
plt.xlim(3000,9500)
plt.ylim(0.,0.8)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)

filtros2 = U.get_str(ruta+'splusNBs.list',0)
plt.subplot(212)
for ii in range(19):
     x,y = U.get_data(filtros2[ii],(0,1))
     plt.plot(x*1000.,y,'-',lw=10,color=colores2[:,ii],alpha=0.5)

for ii in range(19):
     x,y = U.get_data(filtros2[ii],(0,1))
     if ii in NB:
         if ii == 15:
            plt.plot(x*10.,(y/y.max())*0.65,'-',lw=3,color='black',alpha=0.5)
            plt.fill_between(x*10.,(y/y.max())*0.65,0,alpha=0.3,color=colores2[:,ii],lw=4)
         else:
            plt.plot(x*10.,y,'-',lw=3,color='black',alpha=0.5)
            plt.fill_between(x*10.,y,0,alpha=0.3,color=colores2[:,ii],lw=4)
     else:
         plt.plot(x,y*0.65,'-',lw=2,color=colores2[:,ii],alpha=0.5)
         plt.fill_between(x,y*0.65,0,alpha=0.5,color=colores2[:,ii],lw=4)

lab2 = ['J0378','J0395','J0410','J0430','S0446','S0463','S0480','S0497',
          'J0515','S0533','S0553','S0573','S0593','S0613','S0633','J0660',
           'S0676','S0691','J0861']
# plt.legend(lab2,loc='upper right',fontsize=16)
plt.xlabel('Wavelength [$\AA$]',size=28,labelpad=5)
plt.ylabel('Throughput',size=28,labelpad=10)
plt.grid()
plt.xlim(3000,9500)
plt.ylim(0.,0.8)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)


"""