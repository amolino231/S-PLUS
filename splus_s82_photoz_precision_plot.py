__author__ = 'albertomolino'

"""
This routine creates a triple figure showing
the precision of the photo-z in three magnitude bins.

"""

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U
import matplotlib.pyplot as plt
import numpy as N

# Defining a routine

def numberdensity_1(C1,C2):
    """
    As the other ones but fixing the ranges manually.
    USAGE:  matrix,axis2,axis1 = CC_numberdensity_contour(V-B, R-I)
    """
    axis1  = N.arange(0.0,0.5,0.0025) #SPLUS FWHM mosaic
    axis2  = N.arange(0.0,0.5,0.0025)  #SPLUS FWHM mosaic
    matrix = N.zeros((len(axis1),len(axis2)),float)
    for ii in range(len(axis1)-1):
        for jj in range(len(axis2)-1):
            good = N.greater_equal(C1,axis1[ii]) * N.less_equal(C1,axis1[ii+1])
            good *= N.greater_equal(C2,axis2[jj]) * N.less_equal(C2,axis2[jj+1])
            value1,value2 = U.multicompress(good,(C1,C2))
            number = len(value1)
            matrix[ii,jj] = number
    return matrix,axis2,axis1


def numberdensity_2(C1,C2):
    """
    As the other ones but fixing the ranges manually.
    USAGE:  matrix,axis2,axis1 = CC_numberdensity_contour(V-B, R-I)
    """
    axis1  = N.arange(0.0,0.5,0.0025) #SPLUS FWHM mosaic
    axis2  = N.arange(-0.3,0.3,0.0025)  #SPLUS FWHM mosaic
    matrix = N.zeros((len(axis1),len(axis2)),float)
    for ii in range(len(axis1)-1):
        for jj in range(len(axis2)-1):
            good = N.greater_equal(C1,axis1[ii]) * N.less_equal(C1,axis1[ii+1])
            good *= N.greater_equal(C2,axis2[jj]) * N.less_equal(C2,axis2[jj+1])
            value1,value2 = U.multicompress(good,(C1,C2))
            number = len(value1)
            matrix[ii,jj] = number
    return matrix,axis2,axis1



# Root to data
mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/'
root2cats = mainroot + 'splus_cats_NGSL/'
#bpz_cat = root2cats + 'COSMOSeB11new_recal/master.STRIPE82_Photometry.m21_COSMOSeB11new_recal_redu.bpz'
#bpz_cat = root2cats + 'PriorSM/masterBPZ_PriorSM.bpz'
bpz_cat = root2cats + 'COSMOSeB11new_recal/PriorSM/masterBPZ_PriorSM.bpz'


# Reading data
zb,mo,od,zs,chi2,tb = U.get_data(bpz_cat,(1,10,5,9,8,4))

# Doing some cleaning
good = N.less(chi2,30) * N.greater_equal(mo,14)
good *= N.greater_equal(od,0.01)
good *= N.less_equal(zs,0.5)
zb,mo,od,zs,chi2,tb = U.multicompress(good,(zb,mo,od,zs,chi2,tb))

# Defining three samples.
mag_1 = 17
mag_2 = 19
mag_3 = 21
gm1 = N.less_equal(mo,mag_1) #* N.greater(tb,4.5)
gm2 = N.less_equal(mo,mag_2)
gm3 = N.less_equal(mo,mag_3)

# Defining basis for histograms
base_histo = N.arange(-0.2,0.2,0.01)
base_histo_2 = [-0.1,-0.05,0.0,0.05,0.1]
base_histo_3 = [-0.2,-0.1,0.0,0.1,0.2]
base_line = N.arange(0.,1.,0.1)
base_3 = N.arange(0.,0.6,0.1)

# Precision
#sigma = [0.012,0.022,0.030]
#mu = [0.0011,0.0008,0.0002]
#outl = [0.8,1.9,2.9]
#nums = [3644,23070,41470]
sigma = [0.011,0.021,0.028]
mu = [0.0010,0.0013,0.0013]
outl = [0.2,1.2,2.7]

# PLot
plt.figure(1, figsize=(18,10),dpi=70, facecolor='w', edgecolor='k')
plt.clf()

# ===== Left column  =====
##
uno  = plt.axes([.05,.31,.3,.65])
mat1,ax2,ax1 = numberdensity_1(zb[gm1],zs[gm1])
plt.contourf(ax2,ax1,N.log10(mat1.T),50,linewidths=5,vmin=0.0001,vmax=1.5)
#plt.xlabel(r'$z_{spec}$',size=20)
plt.grid()
plt.ylabel(r'$z_{ph}$',size=28,labelpad=-12)
#plt.setp(uno,title='$r<%i$'%(mag_1))
plt.plot(base_line,base_line,'--',color='white',lw=2)
plt.xticks(base_3,['','','','',''],fontsize=15)
plt.ylim(0.01,0.49)
plt.xlim(0.01,0.49)
linea = '$\sigma_{z}=%.3f$, $\mu=%.3f$, $\eta=%.1f$'%(sigma[0],mu[0],outl[0])
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')

##
dos  = plt.axes([.05,0.10,0.3,0.20])
mat1,ax1,ax2 = numberdensity_2(zs[gm1],(zs[gm1]-zb[gm1])/(1.+zs[gm1]))
plt.contour(ax2,ax1,N.log10(mat1.T),500,linewidths=3,vmin=0.0001,vmax=1.5)
plt.plot(base_line,base_line*0.,'--',color='white',lw=2)
plt.ylabel(r'$\delta z$',size=24,labelpad=1)
plt.xlabel(r'$z_{sp}$',size=28,labelpad=15)
plt.grid()
plt.yticks(base_histo_2,['','-0.05','0.00','0.05',''],size=15)
plt.setp(dos,ylim=(-0.09,0.09),xlim=(0.0,0.5))
plt.xticks(base_3,['','0.1','0.2','0.3','0.4'],fontsize=15)
plt.ylim(-0.09,0.09)
plt.xlim(0.01,0.49)

##
tres = plt.axes([0.06,.735,0.1,0.2])
#tres = plt.axes([0.07,.72,0.1,0.2])
aa,bb,cc = plt.hist((zs[gm1]-zb[gm1])/(1.+zs[gm1]),base_histo,facecolor='blue',alpha=0.2)
aa,bb,cc = plt.hist((zs[gm1]-zb[gm1])/(1.+zs[gm1]),base_histo,histtype='step',color='blue',lw=2)
plt.xticks(base_histo_3,['','-0.1','0.0','0.1',''],size=15)
plt.yticks(fontsize=13)
plt.yticks([])
plt.xlabel('$\delta z$',size=22,labelpad=10)
plt.grid()
plt.setp(tres,xlim=(-0.19,0.19))

# ===== Center column =====
##
cuatro = plt.axes([.35,.31,.3,.65])
mat1,ax2,ax1 = numberdensity_1(zb[gm2],zs[gm2])
plt.contourf(ax2,ax1,N.log10(mat1.T),50,linewidths=5,vmin=0.0001,vmax=1.5)
plt.grid()
#plt.ylabel(r'$z_{ph}$',size=28,labelpad=-12)
plt.setp(cuatro,xlim=(0.01,0.39),ylim=(0.01,0.39))
plt.plot(base_line,base_line,'--',color='white',lw=2)
plt.xticks(base_3,['','','','',''],fontsize=15)
plt.yticks(base_3,['','','','',''],fontsize=15)
plt.ylim(0.01,0.49)
linea = '$\sigma_{z}=%.3f$, $\mu=%.3f$, $\eta=%.1f$'%(sigma[1],mu[1],outl[1])
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')

##
cinco  = plt.axes([.35,0.10,0.3,0.20])
mat1,ax1,ax2 = numberdensity_2(zs[gm2],(zs[gm2]-zb[gm2])/(1.+zs[gm2]))
plt.contour(ax2,ax1,N.log10(mat1.T),500,linewidths=3,vmin=0.0001,vmax=1.5)
plt.plot(base_line,base_line*0.,'--',color='white',lw=2)
plt.xlabel(r'$z_{sp}$',size=28,labelpad=15)
plt.xticks(base_3,['','0.1','0.2','0.3','0.4'],fontsize=15)
plt.grid()
plt.yticks(base_histo_2,['','','','',''],size=15)
plt.setp(cinco,ylim=(-0.09,0.09),xlim=(0.0,0.5))
plt.ylim(-0.09,0.09)
plt.xlim(0.01,0.49)

##
seis = plt.axes([0.36,.735,0.1,0.2])
aa,bb,cc = plt.hist((zs[gm2]-zb[gm2])/(1.+zs[gm2]),base_histo,facecolor='blue',alpha=0.2)
aa,bb,cc = plt.hist((zs[gm2]-zb[gm2])/(1.+zs[gm2]),base_histo,histtype='step',color='blue',lw=2)
plt.xticks(base_histo_3,['','-0.1','0.0','0.1',''],size=15)
plt.yticks(fontsize=13)
plt.yticks([])
plt.xlabel('$\delta z$',size=22,labelpad=10)
plt.grid()
plt.setp(seis,xlim=(-0.19,0.19))

# ===== Right column. =====
##
siete = plt.axes([.65,0.31,.3,.65])
mat1,ax2,ax1 = numberdensity_1(zb[gm3],zs[gm3])
plt.contourf(ax2,ax1,N.log10(mat1.T),50,linewidths=5,vmin=0.0001,vmax=1.5)
plt.grid()
#plt.ylabel(r'$z_{ph}$',size=28,labelpad=-12)
plt.plot(base_line,base_line,'--',color='white',lw=2)
plt.xticks(base_3,['','','','',''],fontsize=15)
plt.setp(siete,xlim=(0.01,0.49),ylim=(0.0,0.5))
plt.yticks(base_3,['','','','',''],fontsize=15)
plt.ylim(0.01,0.49)
linea = '$\sigma_{z}=%.3f$, $\mu=%.3f$, $\eta=%.1f$'%(sigma[2],mu[2],outl[2])
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')
plt.annotate(linea,xy=(0.135,0.02),fontsize=23,color='black')


##
ocho = plt.axes([.65,0.10,0.3,0.20])
mat1,ax1,ax2 = numberdensity_2(zs[gm3],(zs[gm3]-zb[gm3])/(1.+zs[gm3]))
plt.contour(ax2,ax1,N.log10(mat1.T),500,linewidths=3,vmin=0.0001,vmax=1.5)
plt.plot(base_line,base_line*0.,'--',color='white',lw=2)
plt.xticks(base_3,['','0.1','0.2','0.3','0.4'],fontsize=15)
plt.grid()
plt.yticks(base_histo_2,['','','','',''],size=15)
plt.xlabel(r'$z_{sp}$',size=28,labelpad=15)
#plt.setp(ocho,ylim=(-0.09,0.09),xlim=(0.0,0.5))
plt.ylim(-0.09,0.09)
plt.xlim(0.01,0.49)

##
nueve = plt.axes([0.66,0.735,0.1,0.2])
aa,bb,cc = plt.hist((zs[gm3]-zb[gm3])/(1.+zs[gm3]),base_histo,facecolor='blue',alpha=0.2)
aa,bb,cc = plt.hist((zs[gm3]-zb[gm3])/(1.+zs[gm3]),base_histo,histtype='step',color='blue',lw=2)
plt.xticks(base_histo_3,['','-0.1','0.0','0.1',''],size=15)
plt.yticks(fontsize=13)
plt.yticks([])
plt.xlabel('$\delta z$',size=22,labelpad=10)
plt.grid()
plt.setp(nueve,xlim=(-0.19,0.19))


#plt.annotate('$\sigma_{z}$',xy=(,),fontsize=22,color='black')
#plt.annotate('$\sigma_{z}$',xy=(,),fontsize=22,color='black')
#plt.annotate('$\sigma_{z}$',xy=(,),fontsize=22,color='black')