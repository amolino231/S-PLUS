__author__ = 'albertomolino'

## This routines serves to create a master file
## from all individual .flux_comparisons.

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import alhambra_photools as A

onlytype = 1

root = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/'
root += 'S82/Dec2017/splus_cats_NGSL/'
"""
# root to master-cat
root_to_mastercat =  root + 'cat/'
if not os.path.exists(root_to_mastercat):
   cmd = '/bin/mkdir %s '%(root_to_mastercat)
   os.system(cmd)

# root to master-bpz
root_to_masterbpz =  root + 'bpz/'
if not os.path.exists(root_to_masterbpz):
   cmd = '/bin/mkdir %s '%(root_to_masterbpz)
   os.system(cmd)

# root to master-flux.
root_to_masterflux =  root + 'flux/'
if not os.path.exists(root_to_masterflux):
   cmd = '/bin/mkdir %s '%(root_to_masterflux)
   os.system(cmd)

# Master cat
master_cat_cat = root_to_mastercat + 'master.STRIPE82_Photometry.m21.cat'
lista_cats = master_cat_cat+'.list'
cmd = 'ls %s*.spz.*.cat > %s'%(root,lista_cats)
os.system(cmd)
print 'Creating master CAT file.'
if not os.path.exists(master_cat_cat):
   A.appendlistcatalog(lista_cats,master_cat_cat)

# Master BPZ
master_bpz_cat = root_to_masterbpz + 'master.STRIPE82_Photometry.m21.bpz'
lista_bpzs = master_bpz_cat+'.list'
cmd = 'ls %s*.bpz > %s'%(root,lista_bpzs)
os.system(cmd)
print 'Creating master BPZ file.'
if not os.path.exists(master_bpz_cat):
   A.appendlistcatalog(lista_bpzs,master_bpz_cat)

# Master FLUX_COMP.
master_flux_cat = root_to_masterflux + 'master.STRIPE82_Photometry.m21.flux_comparison'
lista_fluxes = master_flux_cat+'.list'
cmd = 'ls %s*.flux_comparison > %s'%(root,lista_fluxes)
os.system(cmd)
print 'Creating master BPZ/FLUXCOMP. file.'
if not os.path.exists(master_flux_cat):
    A.appendlistcatalog(lista_fluxes,master_flux_cat)
"""

if onlytype:
   #root_ot = root + '/OT/'
   #root_ot = root + 'OT_SEDrecal_sorted/'
   root_ot = root + 'COSMOSeB11new_OT/'
   # Master FLUX_COMP.
   master_flux_OT_cat = root_ot + 'master.STRIPE82_Photometry.m21.OT.flux_comparison'
   lista_fluxes = master_flux_OT_cat+'.list'
   cmd = 'ls %s*.flux_comparison > %s'%(root_ot,lista_fluxes)
   os.system(cmd)
   print 'Creating master OT/FLUXCOMP. file.'
   if not os.path.exists(master_flux_OT_cat):
      A.appendlistcatalog(lista_fluxes,master_flux_OT_cat)

   master_bpz_OT_cat = root_ot + 'master.STRIPE82_Photometry.m21.OT.bpz'
   lista_bpzs = master_bpz_OT_cat+'.list'
   cmd = 'ls %s*.bpz > %s'%(root_ot,lista_bpzs)
   os.system(cmd)
   print 'Creating master OT/BPZs. file.'
   if not os.path.exists(master_bpz_OT_cat):
      A.appendlistcatalog(lista_bpzs,master_bpz_OT_cat)