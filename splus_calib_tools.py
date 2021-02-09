__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import alhambra_photools as A
import coeio as C
import bpz_tools as B
import useful as U

def appendcatalogs(catalog1,catalog2,catalogOUT):

    """
    The task appends catalogs using only the catalog1's header. Catalog1(withheader)+catalog2(woheader)
    The final (composed) catalog is saved as catalogOUT.
    NEEDLESS TO SAY BOTH CATALOGS HAVE TO HAVE THE SAME FORMAT (ROWS&COLUMNS) !!!
    -----

    """

    print 'Reading file1: ',catalog1
    data1 = C.loaddata(catalog1)      # Loading the whole catalog1 content.
    head1 = C.loadheader(catalog1)    # Loading the original header1.
    print 'Reading file2: ',catalog2
    data2 = C.loaddata(catalog2)      # Loading the whole catalog2 content.
    head2 = C.loadheader(catalog2)    # Loading the original header2.

    outcat = catalogOUT
    print outcat

    try:
       nf1 = N.shape(data1)[0]
       nc1 = N.shape(data1)[1]
    except:
       nf1 = 1
       nc1 = N.shape(data1)[0]

    try:
       nf2 = N.shape(data2)[0]
       nc2 = N.shape(data2)[1]
    except:
       nf2 = 1
       nc2 = N.shape(data2)[0]

    print 'Dimensions catalogue_1: ncols: %i, nraws: %i'%(nf1,nc1)
    print 'Dimensions catalogue_2: ncols: %i, nraws: %i'%(nf2,nc2)

    if nc1 == nc2:
       nf = nf1+nf2
       nc = nc1
       newdata = N.zeros((nf,nc),float)

       for ii in range(nf1):
           if nf1<2: newdata[ii,:] = data1[:]
           else: newdata[ii,:] = data1[ii,:]

       for ii in range(nf2):
           if nf2<2: newdata[ii+nf1,:] = data2[:]
           else: newdata[ii+nf1,:] = data2[ii,:]

       C.savedata(newdata,outcat, dir="",header=head1)     # Saving and creating the new catalog.

    else:
       print 'Different number of rows between catalogs. Impossible to append catalogs !!'


def appendlistcatalog(lista,outfile='None'):

    """
    It appends a list of catalogs via appendcatalogs
    """
    # Declaring some variables.
    list = U.get_str(lista,0)
    temp = len(lista.split('/')[-1])
    root = lista[:-temp]
    print 'Number of catalogs to be appended: %i' %(len(list))
    print 'Starting with the appendage...'

    for jj in range(len(list)-1):
        print 'Appending catalog %i/%i...' %(jj+1,len(list)-1)
        ii = jj+1
        if jj == 0:
           catalog1 = list[jj]
           catalog2 = list[ii]
           finalcatalog = root+'temporal.cat' # trunkcat+'_%i%i.cat' %(ff,ii)
           raimundo = finalcatalog
        else:
           catalog1 = raimundo
           catalog2 = list[ii] # trunkcat+'_%i%i.cat' %(ff,ii)
           finalcatalog = root+'temporal2.cat'

        appendcatalogs(catalog1,catalog2,finalcatalog)

        if os.path.exists(root+'temporal2.cat'):
           cmd = ''
           cmd += '/bin/rm %s' %(raimundo)
           os.system(cmd)
           cmd = '/bin/mv %s %s'%(root+'temporal2.cat',root+'temporal.cat')
           os.system(cmd)
           raimundo = root+'temporal.cat'


    # Saving the final catalog.
    if outfile=='None':
       final = lista[:-((len(lista.split('.')[-1]))+1)]+'_appended.cat'
    else: final = outfile
    cmd = '/bin/mv %s %s'%(root+'temporal.cat',final)
    os.system(cmd)
    print 'A new catalog created as ',final



def rest_frame_wavelength(filter,redshift):
    """
    It estimates the effective filter wavelength
    at the rest-frame (z=0).
    This routine is used in "splus_SEDrecal.py"
    ----
    :param filter:
    :param redshift:
    :return:
    """
    return B.effective_wavelength(filter)/(1.+redshift)


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


def compress_bpz_catalogues(catalogue,sample,outname):
    """
    It selects a subsample of sources from
    an input catalogue.

    :param catalogue:
    :return:
    """
    head1 = C.loadheader(catalogue)
    data1 = C.loaddata(catalogue)
    C.savedata(data1[sample,:],outname,dir="",header=head1)


def replacing_nans_catalogs(catalog,newname):
    """

    vars = []
    evars = []
    data = C.loaddata(catalog)
    mags = data[:,vars]
    emags = data[:,evars]

    """

    vars = N.array([15,18,21,24,27,30,33,36,39,42,45,48,51,54,
                    57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108,111])
    evars = vars[:]+1
    s2n_vars = vars[:]+2

    data = C.loaddata(catalog)      # Loading the whole catalog content.
    head = C.loadheader(catalog)    # Loading the original header.
    mm = data[:,vars]
    em = data[:,evars]
    s2n = data[:,s2n_vars]

    nl = len(mm[:,0])    # nl is the number of detections inside every single band.
    nf = len(mm[0,:])    # nf is the number of bands inside the catalog.
    newmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved.
    errmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved.
    new_s2n = U.zeros((nl,nf),float)

    for jj in range(len(vars)):
        for ii in range(nl):
            if abs(mm[ii,jj]) > 60.:
               newmag[ii,jj] = -99.0
               errmag[ii,jj] =  0.00
               new_s2n[ii,jj] = -1
            elif s2n[ii,jj]<0.00001:
                 new_s2n[ii,jj] = 0.
            else:
               newmag[ii,jj] = mm[ii,jj]
               errmag[ii,jj] = em[ii,jj]
               new_s2n[ii,jj] = s2n[ii,jj]

    # New values of mags error overwrites now the original data.
    data[:,vars]  = newmag[:,U.arange(nf)]
    data[:,evars] = errmag[:,U.arange(nf)]
    data[:,s2n_vars] = new_s2n[:,U.arange(nf)]
    C.savedata(data,newname, dir="",header=head) # Saving & creating a new catalog.




def mixing_2bpz_files(bpz_cali_1,bpz_cali_2):
    """
    It opens 2 BPZ files and select the results for each galaxy
    with the highest Odds values.
    """

    root_new = os.path.dirname(bpz_cali_1)+'/'

    # It opens the files as strings.
    raw = open(bpz_cali_1,'r')
    bpz_data_1 = raw.read()
    bpz_data_1 = bpz_data_1.split('\n')
    raw.close()

    raw = open(bpz_cali_2,'r')
    bpz_data_2 = raw.read()
    bpz_data_2 = bpz_data_2.split('\n')
    raw.close()

    # It estimates the number of commented elements in header.
    n_head = 0
    for sss in range(100):
        if bpz_data_1[sss][0] == "#":
           n_head += 1

    # It reads the Odds column from the BPZ files to compare.
    odds_1 = U.get_data(bpz_cali_1,5)
    odds_2 = U.get_data(bpz_cali_2,5)

    # Number of galaxies in the catalogue.
    n_gal = len(odds_1)
    #print 'ngal: ',n_gal

    # It creates a new file picking the results from the best run.
    mix_bpz_file = root_new + 'best.bpz'
    if not os.path.exists(mix_bpz_file):
       bpz_file = open(mix_bpz_file,'w')
       for iii in range(n_gal):
           odds_values = [odds_1[iii],odds_2[iii]]
           pos = N.argmax(odds_values)
           #print pos
           if pos == 0:
               linea = '%s  \n'%(bpz_data_1[iii+n_head])
           else:
               linea = '%s  \n'%(bpz_data_2[iii+n_head])
           #print linea
           bpz_file.write(linea)
       bpz_file.close()




def mixing_bpz_files(bpz_cali_1,bpz_cali_2,bpz_cali_3):
    """
    It opens 3 BPZ files and select the results for each galaxy
    with the highest Odds values.
    """

    root_new = os.path.dirname(bpz_cali_1)+'/'

    # It opens the files as strings.
    raw = open(bpz_cali_1,'r')
    bpz_data_1 = raw.read()
    bpz_data_1 = bpz_data_1.split('\n')
    raw.close()

    raw = open(bpz_cali_2,'r')
    bpz_data_2 = raw.read()
    bpz_data_2 = bpz_data_2.split('\n')
    raw.close()

    raw = open(bpz_cali_3,'r')
    bpz_data_3 = raw.read()
    bpz_data_3 = bpz_data_3.split('\n')
    raw.close()

    # It estimates the number of commented elements in header.
    n_head = 0
    for sss in range(100):
        if bpz_data_1[sss][0] == "#":
           n_head += 1

    # It reads the Odds column from the BPZ files to compare.
    odds_1 = U.get_data(bpz_cali_1,5)
    odds_2 = U.get_data(bpz_cali_2,5)
    odds_3 = U.get_data(bpz_cali_3,5)

    # Number of galaxies in the catalogue.
    n_gal = len(odds_1)
    #print 'ngal: ',n_gal

    # It creates a new file picking the results from the best run.
    mix_bpz_file = root_new + 'best.bpz'
    if not os.path.exists(mix_bpz_file):
       bpz_file = open(mix_bpz_file,'w')
       for iii in range(n_gal):
           odds_values = [odds_1[iii],odds_2[iii],odds_3[iii]]
           pos = N.argmax(odds_values)
           #print pos
           if pos == 0:
               linea = '%s  \n'%(bpz_data_1[iii+n_head])
           elif pos == 1:
               linea = '%s  \n'%(bpz_data_2[iii+n_head])
           else:
               linea = '%s  \n'%(bpz_data_3[iii+n_head])
           #print linea
           bpz_file.write(linea)
       bpz_file.close()






def replace_photo_uncert(catalog,columns):
    """

    """
    data = C.loaddata(catalog)      # Loading the whole catalog content.
    head = C.loadheader(catalog)    # Loading the original header.
    mm = A.get_magnitudes(catalog,columns)
    em = A.get_errmagnitudes(catalog,columns)
    filters = B.get_filter_list(columns)

    nl = len(mm[:,0])    # nl is the number of detections inside every single band.
    nf = len(mm[0,:])    # nf is the number of bands inside the catalog.
    errmag = U.zeros((nl,nf),float)  # Where the new photo errors will be saved.

    for jj in range(nf):
        maglim = B.get_limitingmagnitude(mm[:,jj],em[:,jj],1.,0.25)
        print 'Limiting Magnitude for filter %s: %.3f'%(filters[jj],maglim)
        for ii in range(nl):
            if mm[ii,jj] == -99.:
                errmag[ii,jj] = 0.00
            elif mm[ii,jj] == 99.:
                errmag[ii,jj] = maglim
            else:
                errmag[ii,jj] = em[ii,jj]

    # New values of mags error overwrites now the original data.
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    data[:,evars] = errmag[:,U.arange(nf)]
    finalcatalog = catalog[:-3]+'upp.cat'
    C.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.


def get_usefulcolumns(columns):

    """
    It extracts the vars,evars,posref,zpe,zpo information
    from a columns file.
    vars & evars: mag & emags positions inside the catalog.
    ==========================================
    USAGE:
    columns = 'cluster.columns'
    vars,evars,posref,zpe,zpo = get_usefulcolumns(columns)
    """
    filt = U.get_str(columns,0)
    nf = 0
    for ii in range(len(filt)):
        if filt[ii][-4:] == '.res': nf += 1
        if filt[ii] == 'M_0':
           posM0 = ii
    print 'Number of filters detected... ', nf

    filtref = int(U.get_str(columns,1)[posM0])-1
    rawvars = U.get_str(columns,1,nf)
    vars = U.zeros(nf)
    evars = U.zeros(nf)
    for jj in range(nf):
        vars[jj] = int(rawvars[jj].split(',')[0])-1   # -1 because of columns's notation
        evars[jj] = int(rawvars[jj].split(',')[1])-1  # -1 because of columns's notation
        if vars[jj] == filtref: posref = int(vars[jj])
    zpe,zpo = U.get_data(columns,(3,4),nf)
    vars = vars.astype(int)
    evars = evars.astype(int)
    return vars,evars,posref,zpe,zpo



def num_tiles(tiles):
    """
    :param tiles: a list with Tile numbers
    :return: a list with no duplications.
    """
    tile_list=[]
    for tt in tiles:
        if tt not in tile_list:
            tile_list.append(tt)

    return tile_list



def get_seeing_from_data(fwhm,mr):
    """
    :param fwhm:
    :param mr:
    :return:
    """
    dm = 0.1
    base = N.arange(0.,10,dm)
    stars = N.greater_equal(mr,14)*N.less_equal(mr,19)
    b1,b2=N.histogram(fwhm[stars],base)
    pos = N.argmax(b1)
    seeing = b2[pos]+(dm/2.)

    return seeing

def get_seeing_from_data_pro(fwhm,mr):
    """
    :param fwhm:
    :param mr:
    :return:
    """
    dm = 0.1
    base = N.arange(0.,7,dm)
    stars = N.greater_equal(mr,14)*N.less_equal(mr,18)
    b1,b2=N.histogram(fwhm[stars],base)
    pos = N.argmax(b1)
    seeing = b2[pos]+(dm/2.)
    redu  = N.less_equal(fwhm,U.mean_robust(fwhm[stars])+U.std_mad(fwhm[stars])*2.)
    redu *= N.greater_equal(fwhm,U.mean_robust(fwhm[stars])-U.std_mad(fwhm[stars])*2.)
    redu *= N.greater_equal(mr,14)*N.less_equal(mr,18)

    return seeing,redu

def get_galaxies_from_data(fwhm,mr):
    """
    :param fwhm:
    :param mr:
    :return:
    """
    ff,ss = get_seeing_from_data_pro(fwhm,mr)
    f_max = max(fwhm[ss])
    redu = N.greater_equal(mr,14) * N.greater_equal(fwhm,f_max)
    return redu


def get_complete_magnitude(rmag):

    good = N.less(rmag,30)
    dm = 0.1
    basem = N.arange(10,30,dm)
    b1,b2=N.histogram(rmag[good],basem)
    pos = N.argmax(b1)
    peak = b2[pos]+(dm/2.)

    return peak


def select_jplus_stars(fwhm,mr):

    seeing = get_seeing_from_data(fwhm,mr)
    fw2=fwhm[mr<19]
    stars = N.less(fwhm,seeing+(U.std_mad(fw2)*3.))

    return stars


def select_jplus_galaxies(fwhm,mr):

    seeing = get_seeing_from_data(fwhm,mr)
    fw2=fwhm[mr<19]
    galaxies = N.greater(fwhm,seeing+(U.std_mad(fw2)*3.))

    return galaxies



def get_NMAD_vs_seeing(bpzfile,mmax):

    bpzs = U.get_str(bpzfile,0)
    nb = len(bpzs)
    valor = N.zeros(nb)
    for ii in range(nb):
        zb,zs,mo = U.get_data(bpzs[ii],(1,9,10))
        dz=(zb-zs)/(1.+zs)
        good = N.less(mo,mmax)
        valor[ii]=U.std_mad(dz[good])

    return valor



def get_master_photoBPZcat(catlist):

    """
    From a list of cats and bpz, it creates
    a master cat+bpz catalogue
    ----
    :param catlist:
    :param bpzlist:
    :return: none
    """

    master_list = os.path.dirname(catlist)+'master.BPZ.list'
    master_name = os.path.dirname(catlist)+'master.BPZ.cat'
    filename = open(master_list,'w')
    # header

    cats = U.get_str(catlist,0)
    bpzs = U.get_str(bpzlist,0)
    nc = len(cats)
    nb = len(bpzs)
    if nc != nb:
       print 'Dimensions mismatch!'
       sys.exit()

    for ii in range(nc):
        cat = cats[ii]
        bpz = cats[ii][:-3]+'bpz'
        out = cats[ii]+'BPZ.cat'
        A.appendColorproBpz(cat,bpz,out)

        if os.path.exists(out):
           filename.write(out+' \n')
    filename.close()

    # Compile the master catalogue.
    A.appendlistcatalog(master_name, master_name)


def select_galaxies_JPLUS_CATBPZcatalog(catalog):

    tile_num,obj_id,ra,dec,fw = U.get_data(catalog,(0,1,2,3,4))
    u_m, u_em, g_m, g_em, r_m, r_em = U.get_data(catalog,(5,6,7,8,9,10))
    i_m, i_em, z_m, z_em, f395_m, f395_em = U.get_data(catalog,(11,12,13,14,15,16))
    f410_m, f410_em, f430_m, f430_em, f515_m, f515_em = U.get_data(catalog,(17,18,19,20,21,22))
    f660_m, f660_em, f861_m, f861_em = U.get_data(catalog,(23,24,25,26))
    zb,zb1,zb2,tb,ods = U.get_data(catalog,(27,28,29,30,31))

    gs = select_jplus_galaxies(fw,r_m)

    out_cat = catalog[:-3]+'gal.cat'
    header_cat = '# 1.tile_id 2.object_id 3.ra 4.dec 5.fwhm 6.uJAVA 7.uJAVA_err 8.gSDSS \
 9.gSDSS_err 10.rSDSS 11.rSDSS_err 12.iSDSS 13.iSDSS_err 14.zSDSS 15.zSDSS_err \
 16.J0395 17.J0395_err 18.J0410 19.J0410_err 20.J0430 21.J0430_err \
 17.J0515 18.J0515_err 19.J0660 20.J0660_err 21.J0861 22.J0861_err \
 23.zb_peak 24.zb_min_1sig 25.zb_max_1sig 26.Spectral-class 27.Odds'

    U.put_data(out_cat,(tile_num[gs],obj_id[gs],ra[gs],dec[gs],fw[gs],
                u_m[gs], u_em[gs], g_m[gs], g_em[gs], r_m[gs], r_em[gs],
                i_m[gs], i_em[gs], z_m[gs], z_em[gs], f395_m[gs], f395_em[gs],
                f410_m[gs], f410_em[gs], f430_m[gs], f430_em[gs],
                f515_m[gs], f515_em[gs], f660_m[gs], f660_em[gs],
                f861_m[gs], f861_em[gs],
                zb[gs],zb1[gs],zb2[gs],tb[gs],ods[gs]),header_cat)




def select_galaxies_JPLUS_catalog(catalog):
    """

import sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
import jplus_calib_tools
from jplus_calib_tools import select_jplus_galaxies,num_tiles
catalog = '/Users/albertomolino/jplus_data_download/SV02_March07/SV02_March07.clean.cat'

    """
    

    tile_num,obj_id,ra,dec,fw = U.get_data(catalog,(0,1,2,3,4))
    u_m, u_em, g_m, g_em, r_m, r_em = U.get_data(catalog,(5,6,7,8,9,10))
    i_m, i_em, z_m, z_em, f395_m, f395_em = U.get_data(catalog,(11,12,13,14,15,16))
    f410_m, f410_em, f430_m, f430_em, f515_m, f515_em = U.get_data(catalog,(17,18,19,20,21,22))
    f660_m, f660_em, f861_m, f861_em = U.get_data(catalog,(23,24,25,26))

    header_cat = '# 1.tile_id 2.object_id 3.ra 4.dec 5.fwhm 6.uJAVA 7.uJAVA_err 8.gSDSS \
    9.gSDSS_err 10.rSDSS 11.rSDSS_err 12.iSDSS 13.iSDSS_err 14.zSDSS 15.zSDSS_err \
    16.J0395 17.J0395_err 18.J0410 19.J0410_err 20.J0430 21.J0430_err \
    22.J0515 23.J0515_err 24.J0660 25.J0660_err 26.J0861 27.J0861_err '
    
    only_tiles = num_tiles(tile_num) # tiles w/o duplications
    n_tiles = len(only_tiles) # number of Tiles

    for ii in range(n_tiles):
        ref_tile = only_tiles[ii]
        g1=N.less(abs(tile_num-ref_tile),1)
        tile_cat = catalog[:-3]+'%i.gal.cat'%(ref_tile)
        tile_num_r,obj_id_r,ra_r,dec_r,fw_r = U.multicompress(g1,(tile_num,obj_id,ra,dec,fw))
        u_m_r, u_em_r, g_m_r, g_em_r, r_m_r, r_em_r = U.multicompress(g1,(u_m, u_em, g_m, g_em, r_m, r_em))
        i_m_r, i_em_r, z_m_r, z_em_r, f395_m_r, f395_em_r = U.multicompress(g1,(i_m, i_em, z_m, z_em, f395_m, f395_em))
        f410_m_r, f410_em_r, f430_m_r, f430_em_r, f515_m_r, f515_em_r = U.multicompress(g1,(f410_m, f410_em, f430_m, f430_em, f515_m, f515_em))
        f660_m_r, f660_em_r, f861_m_r, f861_em_r = U.multicompress(g1,(f660_m, f660_em, f861_m, f861_em))
        gs=select_jplus_galaxies(fw_r,r_m_r)
        # out_cat = maincat[:-3]+'gal.cat'
        U.put_data(tile_cat,(tile_num_r[gs],obj_id_r[gs],ra_r[gs],dec_r[gs],fw_r[gs],u_m_r[gs],
                             u_em_r[gs], g_m_r[gs], g_em_r[gs], r_m_r[gs], r_em_r[gs],i_m_r[gs],
                             i_em_r[gs], z_m_r[gs], z_em_r[gs], f395_m_r[gs], f395_em_r[gs]
                             ,f410_m_r[gs], f410_em_r[gs], f430_m_r[gs], f430_em_r[gs]
                             ,f515_m_r[gs], f515_em_r[gs], f660_m_r[gs], f660_em_r[gs]
                             ,f861_m_r[gs], f861_em_r[gs]),header_cat)




def minimum_photouncert(catalog,columns):
    """

    """
    data = C.loaddata(catalog)      # Loading the whole catalog content.
    head = C.loadheader(catalog)    # Loading the original header.
    mm = A.get_magnitudes(catalog,columns)
    em = A.get_errmagnitudes(catalog,columns)
    nl = len(mm[:,0])    # nl is the number of detections inside every single band.
    nf = len(mm[0,:])    # nf is the number of bands inside the catalog.
    errmag = N.zeros((nl,nf),float)  # Where the new photo errors will be saved.

    for jj in range(nf):
        for ii in range(nl):
            if em[ii,jj] < 0.01: errmag[ii,jj] = 0.03
            elif em[ii,jj] > 1.0: errmag[ii,jj] = 1.0
            else:  errmag[ii,jj] = em[ii,jj]

    # New values of mags error overwrites now the original data.
    vars,evars,posref,zpe,zpo = A.get_usefulcolumns(columns)
    data[:,evars] = errmag[:,N.arange(nf)]
    finalcatalog = catalog[:-3]+'ecor.cat'
    C.savedata(data,finalcatalog, dir="",header=head) # Saving & creating a new catalog.
