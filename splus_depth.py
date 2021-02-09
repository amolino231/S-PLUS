__author__ = 'albertomolino'

# /Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/
import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import matplotlib.pyplot as plt

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/data_quality/'
filename = root+'splus_master.cat'
root_final = root+'Depth_Analysis/'

#Position in the catalogue for each filter
mag_pos = N.array([ 18,  27,  36,  45,  54,  63,  72,  81,  90,  99, 108, 117])
s2n_pos = N.array([ 20,  29,  38,  47,  56,  65,  74,  83,  92,  101, 110, 119])

#colors for different filters
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

# filter names
filters = ['uJAVA','J0378','J0395','J0410','J0430','gSDSS',
          'J0515','rSDSS','J0660','iSDSS','J0861','zSDSS']

# Field-of-View
fov = 1.4*11.

# Base to define the plot
base=N.arange(12,30,.25)
base2 = base[:-1]+((base[1]-base[0])/2.)
base3=N.arange(1,1.0e+6,1)

# Generating figures
for ii in range(12):
    #ii=5
    mr,s2n = U.get_data(filename,(mag_pos[ii],s2n_pos[ii]))
    mr += 0.083 

    c0 = N.greater_equal(s2n,0.)*N.less(mr,30)
    c1 = N.greater_equal(s2n,1.)*N.less(mr,30)
    c2 = N.greater_equal(s2n,3.)*N.less(mr,30)
    c3 = N.greater_equal(s2n,5.)*N.less(mr,30)
    # c4 = N.greater_equal(s2n,100.)*N.less(mr,30)

    plt.figure(12)
    plt.clf()

    a10,a20,a30 = plt.hist(mr[c0],base,alpha=0.5,color='blue')
    a11,a21,a31 = plt.hist(mr[c1],base,alpha=0.5,color='green')
    a12,a22,a32 = plt.hist(mr[c2],base,alpha=0.5,color='orange')
    a13,a23,a33 = plt.hist(mr[c3],base,alpha=0.5,color='red')
    # a14,a24,a34 = plt.hist(mr[c4],base,alpha=0.5,color='red')

    max0=base2[N.argmax(a10)]
    max1=base2[N.argmax(a11)]
    max2=base2[N.argmax(a12)]
    max3=base2[N.argmax(a13)]
    # max4=base2[N.argmax(a14)]

    plt.figure(1)
    plt.subplot(3,4,ii+1)

    plt.plot(base2,a10/fov,'-',lw=8,alpha=0.5,color='blue')
    plt.plot(base2,a11/fov,'-',lw=8,alpha=0.5,color='green')
    plt.plot(base2,a12/fov,'-',lw=8,alpha=0.5,color='orange')
    plt.plot(base2,a13/fov,'-',lw=8,alpha=0.5,color='red')
    # plt.plot(base2,a14/fov,'-',lw=12,alpha=0.5,color='red')
    
    plt.fill_between(base2,a10/fov,0,alpha=0.5,color='blue')
    plt.fill_between(base2,a11/fov,0,alpha=0.5,color='green')
    plt.fill_between(base2,a12/fov,0,alpha=0.5,color='orange')
    plt.fill_between(base2,a13/fov,0,alpha=0.5,color='red')
    # plt.fill_between(base2,a14/fov,0,alpha=0.5,color='red')

    plt.plot(base3*0+max0,base3,'--',lw=1,alpha=0.5,color='blue')
    plt.plot(base3*0+max1,base3,'--',lw=1,alpha=0.5,color='green')
    plt.plot(base3*0+max2,base3,'--',lw=1,alpha=0.5,color='orange')
    plt.plot(base3*0+max3,base3,'--',lw=1,alpha=0.5,color='red')
    # plt.plot(base3*0+max4,base3,'--',lw=3,alpha=0.5,color='red')

    plt.title('%s'%(filters[ii]),size=14)
    plt.xticks(fontsize=15);plt.yticks(fontsize=15)
    if ii>7: plt.xlabel('magnitude',size=17,labelpad=5)
    if ii in [0,4,8]: plt.ylabel('number counts [deg$^{-2}$]',size=16,labelpad=3)
    plt.xlim(13,24)
    plt.ylim(10.,max(a10/fov)*1.2)
    plt.legend(['$S/N$ $>1:$    $%.2f$'%(max0),'$S/N$ $>3:$    $%.2f$'%(max1),'$S/N$ $>5:$    $%.2f$'%(max2),'$S/N$ $>10:$   $%.2f$'%(max3)],loc='upper left',fontsize=14)
    # plt.legend(['$S/N$ $>0:$    $%.1f$'%(max0),'$S/N$ $>3:$    $%.1f$'%(max1),'$S/N$ $>5:$    $%.1f$'%(max2),'$S/N$ $>10:$   $%.1f$'%(max3),'$S/N$ $>100:$ $%.1f$'%(max4),],loc='upper left',fontsize=25)
    plt.grid()
    #pausa = raw_input('paused')
    #plt.savefig(root_final+'dept_%s.png'%(filters[ii]),dpi=80)




