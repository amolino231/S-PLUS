__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import h5py
import matplotlib.pyplot as plt
import tables

def compute_Wittman_CIs_gaussConv(hdf5file,zs,gfactor):
    """
    It calculates the (HPD) CI intervals
    and F(C) for galaxies with zs redshift values.

    """
    verbose=0

    p = h5py.File(hdf5file, mode='r')
    pdz = p.get('FullProbability')
    zz  = p.get('redshift')[:]
    ngal = N.shape(pdz)[0]
    #print 'N.shape(pdz)',N.shape(pdz)
    if ngal <> len(zs):
       print 'Dimension missmatch!'
       sys.exit()
    ci_values = N.zeros(ngal)

    for ii in range(ngal):
        pdz_gal = U.sum(pdz[ii,:,:],axis=1)
        dz = 0.001
        x=N.arange(-3.*gfactor,3.*gfactor+dz/100.,dz)
        # sigma_g=0.02
        gaus = N.exp(-(x/gfactor)**2)
        pdz_gal = N.convolve(pdz_gal,gaus,1)
        pdz_gal /= pdz_gal.sum() #Norm the distr.
        if verbose:
            print 'PDF-size',N.shape(pdz_gal)
            print 'z-length',len(zz)
        dz = abs(zz-zs[ii])
        pos_z = N.where(dz==min(dz))[0][0]
        pdz_th_value = pdz_gal[pos_z]
        if verbose:
           print 'pdz_th_value',pdz_th_value
        good_z_ranges = N.greater_equal(pdz_gal,pdz_th_value)
        pdz_r = N.compress(good_z_ranges,pdz_gal)
        ci_values[ii] = pdz_r.sum()
        #print pdz_r.sum()
        #if ii<10:
        #   plt.clf()
        #   plt.plot(zz,pdz_gal,'k-')
        #   plt.plot(zz[good_z_ranges],pdz_gal[good_z_ranges],'-r.')
        #   plt.plot(zz[pos_z],pdz_gal[pos_z],'go',ms=10)
        #pausa = raw_input('paused')

    return ci_values


def compute_Wittman_CIs(hdf5file,zs):
    """
    It calculates the (HPD) CI intervals
    and F(C) for galaxies with zs redshift values.

    """
    verbose=0

    p = h5py.File(hdf5file, mode='r')
    pdz = p.get('FullProbability')
    zz  = p.get('redshift')[:]
    ngal = N.shape(pdz)[0]
    print 'N.shape(pdz)',N.shape(pdz)
    if ngal <> len(zs):
       print 'Dimension missmatch!'
       sys.exit()
    ci_values = N.zeros(ngal)

    for ii in range(ngal):
        pdz_gal = U.sum(pdz[ii,:,:],axis=1)
        # dz = 0.001
        # x=N.arange(-3.*sigma_g,3.*sigma_g+dz/100.,dz)
        # sigma_g=0.02
        # gaus = N.exp(-(x/sigma_g)**2)
        # pepe = N.convolve(pdz_gal,gaus,1)
        pdz_gal /= pdz_gal.sum() #Norm the distr.
        if verbose:
            print 'PDF-size',N.shape(pdz_gal)
            print 'z-length',len(zz)
        dz = abs(zz-zs[ii])
        pos_z = N.where(dz==min(dz))[0][0]
        pdz_th_value = pdz_gal[pos_z]
        if verbose:
           print 'pdz_th_value',pdz_th_value
        good_z_ranges = N.greater_equal(pdz_gal,pdz_th_value)
        pdz_r = N.compress(good_z_ranges,pdz_gal)
        ci_values[ii] = pdz_r.sum()
        #print pdz_r.sum()
        #if ii<10:
        #   plt.clf()
        #   plt.plot(zz,pdz_gal,'k-')
        #   plt.plot(zz[good_z_ranges],pdz_gal[good_z_ranges],'-r.')
        #   plt.plot(zz[pos_z],pdz_gal[pos_z],'go',ms=10)
        #pausa = raw_input('paused')

    return ci_values


def get_PDZerrDistribution(hdf5file,bpzfile,m_max,o_min):
    """
    It returns the error distribution based on PDZs.
---

import splus_s82_hdf5_tools as to
root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
root += 'S82/Dec2017/splus_cats_NGSL/'
hdf5list = root+'hdf5.list'
bpzlist = root+'bpz/master.STRIPE82_Photometry.m21.bpz.list'
hdf5_files = U.get_str(hdf5list,0)
n_hdf5 = len(hdf5_files)
bpz_files  = U.get_str(bpzlist,0)
n_bpz = len(bpz_files)
for ii in range(n_bpz):
    name = os.path.basename(hdf5_files[ii])
    print name
    try: z,dp,df = to.get_PDZerrDistribution(hdf5_files[ii],bpz_files[ii],20)
    except: print 'Impossible to run on ',name

    """

    plots = 0

    try: ids,zb,zs,mo,tb,odd = U.get_data(bpzfile,(0,1,11,12,4,5))
    except: ids,zb,zs,mo,tb,odd = U.get_data(bpzfile,(0,1,9,10,4,5))
    good = N.less_equal(mo,m_max) * N.greater_equal(odd,o_min)
    ids,zb,zs,mo,tb,odd = U.multicompress(good,(ids,zb,zs,mo,tb,odd))
    ng = len(ids)

    #Readin the PDZs...
    p = h5py.File(hdf5file, mode='r')
    #pdzo = p.get('FullProbability')
    pdz = p.get('Likelihood')
    pdz = pdz[good,:,:]
    zz = p.get('redshift')[:]
    dz = (zz[2]-zz[1])
    basez2 = N.arange(-0.1,0.1,dz)
    basez3 = N.arange(-0.1,0.1,dz*10.)
    basez2b = basez2[:-1]+((basez2[1]-basez2[0])/2.)
    basez3b = basez3[:-1]+((basez3[1]-basez3[0])/2.)
    nz = len(basez2)
    res = 1

    # starting plots if necessary
    #if plots:
    #   plt.figure(12, figsize = (8.5,10.),dpi=80, facecolor='w', edgecolor='k')

    # Computing the z error distr. function
    # based on peak values.
    delta_z_peaks=(zb-zs)/(1.+zs)
    a1,a2 = N.histogram(delta_z_peaks,basez3)

    delta_z_pdzs = N.zeros(nz-1)
    print 'ng:',ng

    for ii in range(ng):
        pdz_mot=U.sum(pdz[ii,:,:],axis=1)
        pdz_mot_peak = pdz_mot/float(max(pdz_mot))
        # To get rid of long tails in PDFs with low probabilities.
        pdz_mot_peak = N.where(pdz_mot_peak<1.0e-4,0.,pdz_mot_peak)
        pdz_mot_norm  = pdz_mot_peak/float(sum(pdz_mot_peak))
        #pdz_mot_norm = N.where(pdz_mot_norm<0.,0.,pdz_mot_norm)
        #pdz_mot_norm  = pdz_mot/float(sum(pdz_mot))
        pdz_mot_norm_resample = U.match_resol(zz-zs[ii],pdz_mot_norm,basez2b)
        pdz_mot_norm_resample = N.where(pdz_mot_norm_resample<0.,0.,pdz_mot_norm_resample)
        delta_z_pdzs += pdz_mot_norm_resample[:]


    # New variables to handle data easily.
    # It scales the normalized PDFs by the ng!
#    norm_dz_peaks = a1/float(sum(a1))
#    norm_dz_pdfs = delta_z_pdzs/float(sum(delta_z_pdzs))
    norm_dz_peaks = a1/float(max(a1))
    norm_dz_pdfs = delta_z_pdzs/float(max(delta_z_pdzs))

    if plots:
       plt.figure(11, figsize = (8.5,10.),dpi=80, facecolor='w', edgecolor='k')
       plt.clf()
       #plt.subplot(212)
       plt.plot(basez3b,norm_dz_peaks,'b-',lw=8,alpha=0.6)
       plt.plot(basez2b,norm_dz_pdfs,'r-',lw=5,alpha=0.9)
       plt.grid()
       plt.xlim(-0.2,0.2)
       plt.ylabel('P(z)',size=20,labelpad=+1)
       plt.legend(['peaks','pdfs'],loc='upper left',fontsize=20)
       plt.xlabel('$\delta_{z}$',size=30)
       final_name = os.path.dirname(hdf5file)+'/PDFs/'+os.path.basename(hdf5file)
       plot_filename = final_name[:-4]+'deltaz.mmax%.2fAB.png'%(m_max)
       #plt.savefig(plot_filename,dpi=80)

    # Saving data into a file.
    #output_filename = hdf5file[:-4]+'deltaz.mmax%.2fAB.mat'%(m_max)
    #U.put_data(output_filename,(basez2b,norm_dz_peaks,norm_dz_pdfs),'z dz_peak dz_PDFs')

    return basez2b,basez3b,norm_dz_peaks,norm_dz_pdfs,ng



def master_PDZerrDistribution(hdf5list,bpzlist,m_max,o_min):

    plots = 1

    hdf5_files = U.get_str(hdf5list,0)
    n_hdf5 = len(hdf5_files)
    bpz_files  = U.get_str(bpzlist,0)
    n_bpz = len(bpz_files)

    if n_bpz != n_hdf5:
       print 'Dimensions mismatch!'
       sys.exit()

    for ii in range(n_hdf5):
        print 'Reading file: ', hdf5_files[ii]
        if ii<1:
           basez,basez2,dz_peak,dz_pdf,ng = get_PDZerrDistribution(hdf5_files[ii],bpz_files[ii],m_max,o_min)
        else:
           basez,basez2,dz_peak_temp,dz_pdf_temp,ng_temp = get_PDZerrDistribution(hdf5_files[ii],bpz_files[ii],m_max,o_min)
           dz_peak += dz_peak_temp
           dz_pdf  += dz_pdf_temp
           ng += ng_temp


    if plots:
       plt.figure(12, figsize = (8.5,10.),dpi=80, facecolor='w', edgecolor='k')
       plt.clf()
       plt.plot(basez2,dz_peak,'b-',lw=12,alpha=0.6)
       plt.plot(basez,dz_pdf,'r-',lw=5,alpha=0.9)
       plt.grid()
       plt.xlim(-0.1,0.1)
       plt.ylabel('P(z)',size=20,labelpad=+1)
       plt.legend(['peaks','pdfs'],loc='upper left',fontsize=20)
       plt.xlabel('$\delta_{z}$',size=30)

    # Saving data into a file.
    #output_filename = hdf5list[:-4]+'master.deltaz.mmax%.2fAB.mat'%(m_max)
    #U.put_data(output_filename,(basez,dz_peak,dz_pdf),'z dz_peak dz_PDFs')

    return basez,basez2,dz_peak,dz_pdf,ng




def get_PDZerrDistribution_byTemplates(hdf5file,bpzfile,m_max):
    """
    It returns the error distribution based on PDZs.
---
import splus_s82_hdf5_tools as to
root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
root += 'S82/Dec2017/splus_cats_NGSL/'
hdf5list = root+'hdf5.list'
bpzlist = root+'bpz/master.STRIPE82_Photometry.m21.bpz.list'
hdf5_files = U.get_str(hdf5list,0)
n_hdf5 = len(hdf5_files)
bpz_files  = U.get_str(bpzlist,0)
n_bpz = len(bpz_files)
for ii in range(n_bpz):
    name = os.path.basename(hdf5_files[ii])
    print name
    try: z,dp,df = to.get_PDZerrDistribution_byTemplates(hdf5_files[ii],bpz_files[ii],19)
    except: print 'Impossible to run on ',name

    """

    plots = 1
    # starting plots if necessary
    if plots:
       plt.figure(12, figsize = (8.5,10.),dpi=80, facecolor='w', edgecolor='k')

    try: ids,zb,zs,mo,tb,odd = U.get_data(bpzfile,(0,1,11,12,4,5))
    except: ids,zb,zs,mo,tb,odd = U.get_data(bpzfile,(0,1,9,10,4,5))
    good = N.less_equal(mo,m_max)
    ids,zb,zs,mo,tb,odd = U.multicompress(good,(ids,zb,zs,mo,tb,odd))
    ng = len(ids)

    #Readin the PDZs...
    p = h5py.File(hdf5file, mode='r')
    #pdzo = p.get('FullProbability')
    pdz = p.get('Likelihood')
    pdz = pdz[good,:,:]
    zz = p.get('redshift')[:]
    dz = (zz[2]-zz[1])*100.
    basez2 = N.arange(-0.2,0.2,dz)
    basez2b = basez2[:-1]+((basez2[1]-basez2[0])/2.)
    nz = len(basez2)
    res = 1

    # Computing the z error distr. function
    # based on peak values.
    delta_z_peaks=(zb-zs)/(1.+zs)
    a1,a2 = N.histogram(delta_z_peaks,basez2)

    delta_z_pdzs = N.zeros(nz-1)
    for ii in range(ng):
        pdz_mot=U.sum(pdz[ii,:,:],axis=1)
        pdz_mot_peak = pdz_mot/float(max(pdz_mot))
        # To get rid of long tails in PDFs with low probabilities.
        pdz_mot_peak = N.where(pdz_mot_peak<1.0e-4,0.,pdz_mot_peak)
        pdz_mot_norm  = pdz_mot_peak/float(sum(pdz_mot_peak))
        pdz_mot_norm = N.where(pdz_mot_norm<0.,0.,pdz_mot_norm)
        #pdz_mot_norm  = pdz_mot/float(sum(pdz_mot))
        pdz_mot_norm_resample = U.match_resol(zz-zs[ii],pdz_mot_norm,basez2b)
        pdz_mot_norm_resample = N.where(pdz_mot_norm_resample<0.,0.,pdz_mot_norm_resample)
        delta_z_pdzs += pdz_mot_norm_resample[:]


        """
        if plots:
           plt.clf()
           plt.subplot(121)
           peak_zb_pos = N.argmax(pdz_mot_norm[::res])
           print zz[peak_zb_pos]
           plt.plot(zz[::res]-zs[ii],pdz_mot_norm[::res],'-',lw=5,alpha=0.6)
           #plt.plot(zz[::res]-zz[peak_zb_pos],pdz_mot_norm[::res],'-',lw=5,alpha=0.6)
           plt.grid()
           plt.xlim(-0.2,0.2)
           #plt.ylim(0.001,0.1)
           plt.xlabel('$\delta_{z}$',size=30)
           plt.ylabel('P(z)',size=20,labelpad=+1)
           plt.legend(['R=%.2f''\n''T=%.1f''\n''O=%.1f'%(mo[ii],tb[ii],odd[ii])],loc='upper right')
           plt.title('zb = %.2f, zs = %.2f, dz/1+z = %.2f'%(zb[ii],zs[ii],delta_z_peaks[ii]),size=20)
           plt.subplot(122)
           plt.plot(basez2b,delta_z_pdzs,'k-',lw=5)
           plt.grid()
           plt.xlim(-0.2,0.2)
           #plt.ylim(0.001,0.1)
           plt.xlabel('$\delta_{z}$',size=30)
           plt.ylabel('P(z)',size=20,labelpad=+1)
           pausa = raw_input('press a bottom to continue')
        """


    # New variables to handle data easily.
    # It scales the normalized PDFs by the ng!
    norm_dz_peaks = a1/float(sum(a1))
    norm_dz_pdfs = delta_z_pdzs/float(sum(delta_z_pdzs))

    if plots:
       plt.figure(11, figsize = (8.5,10.),dpi=80, facecolor='w', edgecolor='k')
       plt.clf()
       #plt.subplot(212)
       plt.plot(basez2b,norm_dz_peaks,'b-',lw=8,alpha=0.6)
       plt.plot(basez2b,norm_dz_pdfs,'r-',lw=5,alpha=0.9)
       plt.grid()
       plt.xlim(-0.2,0.2)
       plt.ylabel('P(z)',size=20,labelpad=+1)
       plt.legend(['peaks','pdfs'],loc='upper left',fontsize=20)
       plt.xlabel('$\delta_{z}$',size=30)
       plot_filename = hdf5file[:-4]+'deltaz.mmax%.2fAB.png'%(m_max)
       plt.savefig(plot_filename,dpi=80)

    # Saving data into a file.
    output_filename = hdf5file[:-4]+'deltaz.mmax%.2fAB.mat'%(m_max)
    U.put_data(output_filename,(basez2b,norm_dz_peaks,norm_dz_pdfs),'z dz_peak dz_PDFs')

    return basez2b,norm_dz_peaks,norm_dz_pdfs


def getPDF_by_mag_templates_and_weights(hdf5file,m_max,weights):
    """
    It returns the global P(z) for each template individually.
    ---
    m_mag: maximum magnitude to be considered.
    weights: vector with P(det|gal).

    """
    p = h5py.File(hdf5file, mode='r')
    # pdz = p.get('FullProbability')
    pdz = p.get('Likelihood')
    z   = p.get('redshift')
    mmm = p.get('m_0')[:]
    #good_m_sample = N.less_equal(mmm,m_max)
    #pdz = pdz[good_m_sample,:,:]
    ngal = N.shape(pdz)[0]
    #weights_redu = weights[good_m_sample]
    weights_redu = weights * 1.
    tts  = p.get('type')[:]
    tbs = N.unique(tts.astype(int)) # indiv integers.
    ntb = len(tbs)

    global_pdfs = N.zeros((len(z),ntb),'float')
    ngal = 20
    for ii in range(ngal):
        # I need to separate the PDF per each type.
        print 'Analyzing galaxy %i'%(ii+1)
        for uu in range(ntb):
            if uu<1:
                good_T_sample = N.less_equal(tts,tbs[uu+1]-0.5)
            elif ii == ntb-1:
                good_T_sample = N.greater_equal(tts,tbs[uu]-0.5)
            else:
                good_T_sample = N.greater_equal(tts,tbs[uu]-0.5)
                good_T_sample *= N.less_equal(tts,tbs[uu]+0.5)

            # Global PDF_ii
            global_pdf_galaxy_ii = U.sum(pdz[ii,:,:],axis=1)

            # Now PDF_ii for each template (uu: good_T_sample).
            if ii < 1:
               global_pdfs[:,uu]  = U.sum(pdz[ii,:,good_T_sample],axis=1)
               if global_pdfs[:,uu].sum() > 0.:
                  global_pdfs[:,uu] /= (1.*global_pdf_galaxy_ii) # Normalize
                  global_pdfs[:,uu] *= weights_redu[ii] # Probability of being galaxy.
               else:
                  global_pdfs[:,uu] = N.zeros(len(z))

            else:
               global_pdfs[:,uu] += U.sum(pdz[ii,:,good_T_sample],axis=1)
               if global_pdfs[:,uu].sum() > 0.:
                   global_pdfs[:,uu] /= (1.*global_pdf_galaxy_ii) # Normalize
                   global_pdfs[:,uu] *= weights_redu[ii] # Probability of being galaxy.
               else:
                   global_pdfs[:,uu] += N.zeros(len(z))



    return z,global_pdfs




def getPDF_by_mag_and_weights(hdf5file,m_max,weights):
    """
    It returns the global P(z)
    ---
    m_mag: maximum magnitude to be considered.
    weights: vector with P(det|gal).

    """
    p = h5py.File(hdf5file, mode='r')
    pdz = p.get('FullProbability')
    # pdz = p.get('Likelihood')
    z   = p.get('redshift')
    mmm = p.get('m_0')[:]
    good_sample = N.less_equal(mmm,m_max)
    pdz = pdz[good_sample,:,:]
    ngal = N.shape(pdz)[0]
    weights_redu = weights[good_sample]

    globalpdz_red  = N.zeros(len(z))
    globalpdz_blue = N.zeros(len(z))
    globalpdz3_all = N.zeros(len(z))

    for ii in range(ngal):
        try: pepe1  = U.sum(pdz[ii,:,0:75],axis=1)
        except: pepe1 = z*0.
        try: pepe2  = U.sum(pdz[ii,:,75:],axis=1)
        except: pepe2 = z*0.
        try: pepe66 = U.sum(pdz[ii,:,:],axis=1)
        except: pepe66 = z*0.
        try: pepe3  = (U.sum(pepe1)+U.sum(pepe2))*1.
        except: pepe3 = z*0.

        if ii==0:
           try: globalpdz_red  = (pepe1/pepe3) * weights_redu[ii]
           except: globalpdz_red = z*0.
           try: globalpdz_blue = (pepe2/pepe3) * weights_redu[ii]
           except: globalpdz_blue = z*0.
           # try: globalpdz3_all = pepe66
           try: globalpdz3_all = (pepe66/pepe66.sum()) * weights_redu[ii]
           except: globalpdz_all = z*0.
        else:
           try: globalpdz_red  += (pepe1/pepe3) * weights_redu[ii]
           except: globalpdz_red += z*0.
           try: globalpdz_blue += (pepe2/pepe3) * weights_redu[ii]
           except: globalpdz_blue += z*0.
           try: globalpdz3_all += (pepe66/pepe66.sum()) * weights_redu[ii]
           # try: globalpdz3_all += pepe66
           except: globalpdz3_all += z*0.



    return z,globalpdz_red,globalpdz_blue,globalpdz3_all,ngal




def global_PDZ(hdf5file,m_max):
    """
    It returns the global P(z)
    """
    p = h5py.File(hdf5file, mode='r')
    # pdz = p.get('FullProbability')
    pdz = p.get('Likelihood')
    z   = p.get('redshift')
    mmm = p.get('m_0')[:]
    good_sample = N.less_equal(abs(mmm),m_max)
    pdz = pdz[good_sample,:,:]
    ngal = N.shape(pdz)[0]

    globalpdz_red = N.zeros(len(z))
    globalpdz_blue = N.zeros(len(z))
    globalpdz3_all = N.zeros(len(z))
    for ii in range(ngal):
        try: pepe1  = U.sum(pdz[ii,:,0:75],axis=1)
        except: pepe1 = z*0.
        try: pepe2  = U.sum(pdz[ii,:,75:],axis=1)
        except: pepe2 = z*0.
        try: pepe66 = U.sum(pdz[ii,:,:],axis=1)
        except: pepe66 = z*0.
        try: pepe3  = (U.sum(pepe1)+U.sum(pepe2))*1.
        except: pepe3 = z*0.
        if ii==0:
           try: globalpdz_red  = (pepe1/pepe3)
           except: globalpdz_red = z*0.
           try: globalpdz_blue = (pepe2/pepe3)
           except: globalpdz_blue = z*0.
           # try: globalpdz3_all = pepe66
           try: globalpdz3_all = (pepe66/pepe66.sum())
           except: globalpdz_all = z*0.
        else:
           try: globalpdz_red  += (pepe1/pepe3)
           except: globalpdz_red += z*0.
           try: globalpdz_blue += (pepe2/pepe3)
           except: globalpdz_blue += z*0.
           try: globalpdz3_all += (pepe66/pepe66.sum())
           # try: globalpdz3_all += pepe66
           except: globalpdz_all += z*0.

    return z,globalpdz_red/(ngal*1.),globalpdz_blue/(ngal*1.),globalpdz3_all/(ngal*1.)


def master_global_PDF(hdf5list,m_max):
    """
    This routine serves to extract the final P(z)
    from a list of HDF5 files. A magnitude-cut (m<19)
    is applied.


    :param hdf5list: list of HDF5 files
    :return: zz, p_r,p_b_p_a
    """

    plots=1

    hdf5files = U.get_str(hdf5list,0)
    n_fields = len(hdf5files)
    for ii in range(n_fields):
        if ii<1:
           zz, p_red,p_blue,p_global = global_PDZ(hdf5files[ii],m_max)
        else:
           zz,p_red_temp,p_blue_temp,p_global_temp = global_PDZ(hdf5files[ii],m_max)
           p_red += p_red_temp
           p_blue += p_blue_temp
           p_global += p_global_temp

    if plots:
       plt.figure(12, figsize = (8.5,10.),dpi=80, facecolor='w', edgecolor='k')
       plt.clf()
       plt.plot(zz,p_red,'r-',lw=5,alpha=0.7)
       plt.plot(zz,p_blue,'b-',lw=5,alpha=0.7)
       plt.plot(zz,p_global,'k--',lw=5,alpha=0.7)
       plt.grid()
       plt.xlim(0.,zz.max())
       plt.grid()
       plt.ylabel('P(z)',size=20,labelpad=+1)
       plt.legend(['early','late','all'],loc='upper right',fontsize=20)
       plt.xlabel('$z$',size=30)

    output_filename = hdf5list[:-4]+'master.PDF.mmax%.2fAB.mat'%(m_max)
    U.put_data(output_filename,(zz,p_red,p_blue,p_global),'z P_r P_b P_a')

    return zz,p_red,p_blue,p_global




def get_masterPDZ_pro(hdf5file,mmin,mmax,dm):
    """
    TAKEN DIRECTLY FROM CLASH_TOOLS.py

    This routine derives the master PDZ for
    a sample of magnitudes based on an empirical
    BPZ-HDF5 catalogue.
    ---
    This new version allows the user to change
    the magnitude range from outside.
    ===
import clash_tools as CT
hdf5file = root+'alhambra.spz.hdf5'
mpdz,z,sigz,meanz = CT.get_masterPDZ(hdf5file)

    """
    # Reading data
    p1 = h5py.File(hdf5file, mode='r')
    pdz1 = p1.get('FullProbability')
    zz1  = p1.get('redshift')[:]
    tt1  = p1.get('type')[:]
    mo1  = p1.get('m_0')[:]

    # Defining resolution and other variables.
    basem = N.arange(mmin,mmax+dm,dm)
    nm = len(basem)
    res = 2
    zz1r = zz1[::res]
    dz2  = (zz1r.max()-zz1r.min())
    basez = N.linspace(-dz2,dz2,len(zz1r)*2)
    basez2 = basez+((basez[1]-basez[0])/2.)
    masterpdz = N.zeros((nm,len(basez)),float)
    sigma_pdz = N.zeros(nm)
    meanz = N.zeros(nm)

    # Global P(z)
    weirdpeaks = []
    print 'Number of magnitude-bins: ',nm
    for jj in range(nm):
        # Selecting galaxies within that magnitude bin.
        if jj==0:
           good = N.less_equal(mo1,basem[jj])
           pdz1r = pdz1[good,:,:]
        elif jj==nm-1:
           good = N.greater_equal(mo1,basem[jj-1])
           good *= N.less_equal(mo1,basem[jj])
           pdz1r = pdz1[good,:,:]
        else:
           good =  N.greater_equal(mo1,basem[jj-1])
           good *= N.less_equal(mo1,basem[jj])
           pdz1r = pdz1[good,:,:]

        ng1 = N.shape(pdz1r)[0]
        peaks = []
        print '%i galaxies in magnitude-bin %i '%(ng1,jj+1)
        mo1r = mo1[good]
        for ii in range(ng1):
            # print 'galaxy number %i out of %i '%(ii+1,ng1)
            pdz_ind = N.sum(pdz1r[ii,:,:],axis=1)
            pdz_ind_r = pdz_ind[::res]
            peak = pdz_ind_r.max()
            pos  = N.where(pdz_ind_r==peak)[0][0]
            zpeak = zz1r[pos]
            # new_zmin = len(zz1r)-zpeak
            new_zmin = len(zz1r)-pos-1 # len(zz1r[0:pos])
            min_val = min(pdz_ind_r)
            max_val = max(pdz_ind_r)
            if abs(min_val-max_val)<1.0e-8: pdz_ind_r[:]=zz1r*0
            if zpeak<0.0005: pdz_ind_r[:]=zz1r*0
            # try:
            if zpeak>0.0005:
                masterpdz[jj,new_zmin:new_zmin+len(zz1r)] += pdz_ind_r[:]
            # except: weirdpeaks.append(zpeak)
            if zpeak>0.0005:
                peaks.append(zpeak)

            # plt.figure(200)
            # plt.clf()
            # plt.plot(zz1r,pdz_ind_r,'k-')
            # plt.xlim(0.,0.2)
            # plt.grid()
            # print 'AB,dz/1+z,dz: %.2f,%.3f,%.2f:'%(mo1r[ii],(zpeak-0.044)/1.044,zpeak-0.044)
            # pausa = raw_input('paused')


        masterpdz[jj,:] /= float(masterpdz[jj,:].sum())
        cumasterpdz = U.add.accumulate(masterpdz[jj,:])
        cumasterpdz /= cumasterpdz.max()
        zmin_err_e = U.match_resol(cumasterpdz,basez2,0.17)
        zmax_err_e = U.match_resol(cumasterpdz,basez2,0.83)
        sigma_pdz[jj] = (zmax_err_e-zmin_err_e)
        peak_position = N.where(masterpdz[jj,:]==masterpdz[jj,:].max())[0][0]
        # print 'peak_position',peak_position
        # print 'zpeak=',zz1r[peak_position]
        meanz[jj] = U.mean_robust(N.array(peaks))
        # U.std_mad(N.array(peaks))

    return masterpdz,basez2,sigma_pdz,meanz