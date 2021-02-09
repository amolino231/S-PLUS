__author__ = 'albertomolino'

#/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/

"""
- The motivation of this code is to pick the best photo-z solution
for a given galaxy from different estimates using different photometries.
- The idea is to make the selection based on the photometry that provides
a higher Odds value.
- This tools creates a new BPZ catalogue with the best solutions.

"""

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import useful as U
import numpy as N
#import bpz_tools as B

mainroot = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/S82/Sept2017/'
root2cats = mainroot + 'splus_cats/'
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root2codes = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/'
lista_cats = 'photometry.list'
cats_names = U.get_str(root2cats+lista_cats,0)
n_cats = len(cats_names)
apertures = ['auto','petro','aper']
n_apers = len(apertures)
# splus_columns = root2cats + 'splus11_auto.columns'
sdss_s82_spz_cat = root2cats + 'S82_SDSS.cat'
# spectra = 'eB11.list'
spectra = 'GOSMOSeB11.list'

for ggg in range(n_cats):
    bpz_cali_1  = cats_names[ggg][:-3] + 'spz.%s_cali.bpz'%(apertures[0])
    bpz_cali_2  = cats_names[ggg][:-3] + 'spz.%s_cali.bpz'%(apertures[1])
    bpz_cali_3  = cats_names[ggg][:-3] + 'spz.%s_cali.bpz'%(apertures[2])

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

    # It creates a new file picking the results from the best run.
    mix_bpz_file = cats_names[ggg][:-3] + 'spz.best_cali.bpz'
    if not os.path.exists(mix_bpz_file):
       bpz_file = open(mix_bpz_file,'w')
       for iii in range(n_gal):
           odds_values = [odds_1[iii],odds_2[iii],odds_3[iii]]
           pos = N.argmax(odds_values)
           if pos == 0:
           elif pos == 1:
           else:
           linea = '%s  \n'%(bpz_data_1[iii+n_head])
           bpz_file.write(linea)
    bpz_file.close()



def mixing_bpz_files(bpz_cali_1,bpz_cali_2,bpz_cali_3):


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
    print 'ngal: ',n_gal

    # It creates a new file picking the results from the best run.
    mix_bpz_file = root_new + 'best.bpz'
    if not os.path.exists(mix_bpz_file):
       bpz_file = open(mix_bpz_file,'w')
       for iii in range(n_gal):
           odds_values = [odds_1[iii],odds_2[iii],odds_3[iii]]
           pos = N.argmax(odds_values)
           print pos
           if pos == 0:
               linea = '%s  \n'%(bpz_data_1[iii+n_head])
           elif pos == 1:
               linea = '%s  \n'%(bpz_data_2[iii+n_head])
           else:
               linea = '%s  \n'%(bpz_data_3[iii+n_head])
           print linea
           bpz_file.write(linea)
       bpz_file.close()


