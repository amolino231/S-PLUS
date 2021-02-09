__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import matplotlib.pyplot as plt

root_to_bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root_to_filts = root_to_bpz+'FILTER/'
filters = U.get_str(root_to_filts+'splusNBs.list',0)

narrow = []
narrow.append('F378')
narrow.append('F395')
narrow.append('F410')
narrow.append('F430')
narrow.append('F515')
narrow.append('F660')
narrow.append('F861')

#colors for different filters
colores = N.zeros((3,7),float)
colores[:,0]=(0.00,0.25,1.00)
colores[:,1]=(0.00,0.65,1.00)
colores[:,2]=(0.00,0.50,0.00)
colores[:,3]=(0.85,0.65,0.00)
#colores[:,5]=(0.75,0.50,0.00)
colores[:,4]=(0.80,0.25,0.00)
#colores[:,7]=(1.00,0.00,0.00)
colores[:,5]=(0.85,0.00,0.00)
#colores[:,9]=(0.65,0.00,0.00)
colores[:,6]=(0.35,0.00,0.00)
#colores[:,11]=(0.25,0.00,0.00)

# Emission lines.
galw = N.zeros(7)
galw[0]=3727
galw[1]=4862
galw[2]=5007
galw[3]=6300
galw[4]=6549
galw[5]=6600
galw[6]=6717
galine = []
galine.append('[OII]: 3727$\AA$')
galine.append('H-beta: 4862$\AA$')
galine.append('[OIII]: 5007$\AA$')
galine.append('[OI]: 6300$\AA$')
galine.append('[NII]: 6549$\AA$')
galine.append('H-alpha: 6600$\AA$')
galine.append('[SII]: 6717$\AA$')
#QSOs
galq = N.zeros(9)
galq[0]=1215
galq[1]=1240
galq[2]=1549
galq[3]=1908
galq[4]=2799
galq[5]=4341
galq[6]=4862
galq[7]=5000
galq[8]=6564

# Redshift
z = N.arange(0.,1.,0.01)

plt.figure(1, figsize=(11,9),dpi=80, facecolor='w', edgecolor='k')
plt.clf()
for ii in range(7):
        x,y = U.get_data(root_to_filts+filters[ii],(0,1))
        plt.plot(x,y,'-',lw=1,color=colores[:,ii])
        plt.fill_betweenx(y,x,0,color=colores[:,ii],alpha=0.1)
        plt.grid()
        plt.xlabel('Wavelength [$\AA$]',size=25,labelpad=5)
        plt.ylabel('$z$',size=35)
        plt.xlim(3500,9500);plt.ylim(0.01,1.)
plt.plot(galw[2]*(1+z),z,'--',color='black',alpha=0.95)
plt.plot(galw[5]*(1+z),z,'-',color='black',alpha=0.95)
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)



"""
qsoline = []
qsoline.append('Ly-alpha: 1215$\AA$')
qsoline.append('$N_{V}$: 1240$\AA$')
qsoline.append('$C_{IV}$: 1549$\AA$')
qsoline.append('$C_{III}$: 1908$\AA$')
qsoline.append('$Mg_{II}$: 2799$\AA$')
qsoline.append('H-gamma: 4341$\AA$')
qsoline.append('H-beta: 4862$\AA$')
qsoline.append('[OIII]: 5000$\AA$')
qsoline.append('H-alpha: 6564$\AA$')


for hh in range(9):
    plt.figure(10, figsize=(11,9),dpi=80, facecolor='w', edgecolor='k')
    plt.clf()
    for ii in range(7):
        x,y = U.get_data(root2filts+'/'+narrow[ii],(0,1))
        plt.plot(x,y/y.max() * 7,'-',lw=1,color=matcol2[:,ii])
        plt.fill_betweenx(y/y.max() * 7,x,0,color=matcol2[:,ii],alpha=0.1)
        plt.grid()
        plt.xlabel('Wavelength [$\AA$]',size=20)
        plt.ylabel('$z$',size=30)
        plt.xlim(3500,9500);plt.ylim(0.01,7.)
    fff = raw_input('paused')
    plt.plot(galq[hh]*(1+z),z,'k--')
    plt.title(qsoline[hh])
    plt.savefig('/Users/albertomolino/codigos/bpz-1.99.2/FILTER/JPLUS_FS_20130118/ELQd.%i.png'%(galq[hh]),dpi=150)



"""