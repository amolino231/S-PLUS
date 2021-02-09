__author__ = 'albertomolino'

# /Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/
import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U
import matplotlib.pyplot as plt

root = '/Users/albertomolino/Desktop/SPLUS_cats/'
filename = root+'master_Stripe82.cat'
root_final = root+'Depth_Analysis/'

#Position in the catalogue for each filter
mag_pos = N.array([ 18,  27,  36,  45,  54,  63,  72,  81,  90,  99, 108, 117])
s2n_pos = N.array([ 20,  29,  38,  47,  56,  65,  74,  83,  92,  101, 110, 119])

mag_3arc_pos = N.array([ 21,  30,  39,  48,  57,  66,  75,  84,  93,  102, 111, 120])
s2n_3arc_pos = N.array([ 23,  32,  41,  50,  59,  68,  77,  86,  95,  104, 113, 122])

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
base  = N.arange(12,25,.25)
base2 = base[:-1]+((base[1]-base[0])/2.)
base3 = N.arange(1,1.0e+6,1)

# Generating figures
for ii in range(12):
    mr,s2n = U.get_data(filename,(mag_pos[ii],s2n_pos[ii]))
    mr += 0.083 
    c3 = N.greater_equal(s2n,3.)*N.less(mr,30)

    plt.figure(1)
    plt.clf()
    plt.subplot(211)
    a13_p,a23_p,a33_p = plt.hist(mr[c3],base,alpha=0.5,color='red')
    m_peak = base2[N.argmax(a13_p)]

    plt.figure(2)
    plt.clf()
    c0 = N.greater_equal(s2n,.01)*N.less(mr,30)
    a13_c,a23_c,a33_c = plt.hist(mr[c0],base,cumulative='yes')
    aa_norm_c = (a13_c/fov)/max(a13_c/fov)
    plt.figure(1)
    plt.subplot(212)
    plt.plot(base2,aa_norm_c,'o-',lw=12,alpha=0.5,color='orange')
    pos_50 = N.argmin(abs(aa_norm_c-0.5))
    m_50 = base2[pos_50]
    pos_80 = N.argmin(abs(aa_norm_c-0.8))
    m_80 = base2[pos_80]
    pos_95 = N.argmin(abs(aa_norm_c-0.95))
    m_95 = base2[pos_95]

    m_c3,s2n_c3 = U.get_data(filename,(mag_3arc_pos[ii],s2n_3arc_pos[ii]))
    m_c3 += -0.2
    c3_min = N.greater_equal(s2n_c3,2.95)*N.less_equal(s2n_c3,3.05)
    c3_min *= N.less(m_c3,30)
    m_3s_3arcs = U.mean_robust(m_c3)

    print '%s &  %.2f  &  %.2f  &  %.2f  & %.2f  & %.2f '%(filters[ii],m_peak,m_50,m_80,m_95,m_3s_3arcs)
    
    pausa = raw_input('paused')
