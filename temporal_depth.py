__author__ = 'albertomolino'

# /Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/
import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import matplotlib.pyplot as plt

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Dec2017/splus_cats/'
root_final = root+'Depth_Analysis/'
if not os.path.exists(root_final):
   cmd = '/bin/mkdir %s '%(root_final)
   os.system(cmd)

# Reading catalogues
lista_catalogues = root+'photometry.list'
catalogues = U.get_str(lista_catalogues,0)
n_cats = len(catalogues)
print 'analyzing %i catalogues '%(n_cats)

#Position in the catalogue for each filter
mag_pos = N.array([ 18,  27,  36,  45,  54,  63,  72,  81,  90,  99, 108, 117])
s2n_pos = N.array([ 20,  29,  38,  47,  56,  65,  74,  83,  92,  101, 110, 119])

# filter names
filters = ['uJAVA','J0378','J0395','J0410','J0430','gSDSS',
          'J0515','rSDSS','J0660','iSDSS','J0861','zSDSS']

# Field-of-View
fov = 1.4*n_cats

# Base to define the plot
base=N.arange(12,30,.25)
base2 = base[:-1]+((base[1]-base[0])/2.)
base3=N.arange(1,1.0e+6,1)

final_a10 = N.zeros(len(base2))
final_a11 = N.zeros(len(base2))
final_a12 = N.zeros(len(base2))
final_a13 = N.zeros(len(base2))

for ii in range(12):
    print 'Analyzing filter: ',filters[ii]
    for hhh in range(n_cats):
        mr,s2n = U.get_data(catalogues[hhh],(mag_pos[ii],s2n_pos[ii]))
        mr += 0.083

        c0 = N.greater_equal(s2n,3.)*N.less(mr,30)
        c1 = N.greater_equal(s2n,5.)*N.less(mr,30)
        c2 = N.greater_equal(s2n,10.)*N.less(mr,30)
        c3 = N.greater_equal(s2n,50.)*N.less(mr,30)

        plt.figure(12)
        plt.clf()
        a10,a20,a30 = plt.hist(mr[c0],base,alpha=0.5,color='blue')
        a11,a21,a31 = plt.hist(mr[c1],base,alpha=0.5,color='green')
        a12,a22,a32 = plt.hist(mr[c2],base,alpha=0.5,color='orange')
        a13,a23,a33 = plt.hist(mr[c3],base,alpha=0.5,color='red')

        final_a10[:] += a10[:]
        final_a11[:] += a11[:]
        final_a12[:] += a12[:]
        final_a13[:] += a13[:]

    max0=base2[N.argmax(final_a10)]
    max1=base2[N.argmax(final_a11)]
    max2=base2[N.argmax(final_a12)]
    max3=base2[N.argmax(final_a13)]

    plt.figure(1)
    plt.subplot(3,4,ii+1)

    plt.plot(base2,final_a10/fov,'-',lw=8,alpha=0.5,color='blue')
    plt.plot(base2,final_a11/fov,'-',lw=8,alpha=0.5,color='green')
    plt.plot(base2,final_a12/fov,'-',lw=8,alpha=0.5,color='orange')
    plt.plot(base2,final_a13/fov,'-',lw=8,alpha=0.5,color='red')

    plt.fill_between(base2,final_a10/fov,0,alpha=0.5,color='blue')
    plt.fill_between(base2,final_a11/fov,0,alpha=0.5,color='green')
    plt.fill_between(base2,final_a12/fov,0,alpha=0.5,color='orange')
    plt.fill_between(base2,final_a13/fov,0,alpha=0.5,color='red')

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
    plt.ylim(10.,2.0e+5)
    plt.legend(['$S/N>3:$ $%.2f$'%(max0),'$S/N>5:$ $%.2f$'%(max1),
                '$S/N>10:$ $%.2f$'%(max2),'$S/N>50:$ $%.2f$'%(max3)]
               ,loc='lower left',fontsize=12)
    plt.grid()
    #pausa = raw_input('paused')
    #plt.savefig(root_final+'dept_%s.png'%(filters[ii]),dpi=80)


"""
for ii in range(12):
    plt.subplot(3,4,ii+1)
    if ii<4: plt.ylim(1.,2.0e+3)
    elif ii>3 and ii<8: plt.ylim(1.,3.0e+3)
    else: plt.ylim(1.,2.0e+4)

"""

