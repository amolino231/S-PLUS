__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
sys.path.append('/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/')
import useful as U
import bpz_tools as B

# General roots.
root = '/Volumes/CLASH/S82/specz/'
root2 = root + 'AbsMag/'
root2bpz = '/Users/albertomolino/codigos/bpz-1.99.2/'
root2codes = '/Users/albertomolino/Postdoc/T80S_Pipeline/Commisioning/codes/'
# Files
catalog = root + 'Master.cat'
columns = root + 'Master.columns'

# Filters
filters = B.get_filter_list(columns)
nf = len(filters)

for ii in range(nf):
    ref_filter = filters[ii]
    print 'Processing filter: ',ref_filter
    filter_nick_name = ref_filter[6:][:-4]
    final_folder = root2 + '%s/'%(filter_nick_name)
    if not os.path.exists(final_folder):
       cmd2 = '/bin/mkdir %s '%(final_folder)
       os.system(cmd2)

    #Remove temporal files from BPZ
    rm_ff = root2codes+'*COSMOSeB11new_recal*'
    print rm_ff
    #pausa = raw_input('paused')
    try:
       cmd1 = '/bin/rm %s '%(rm_ff)
       os.system(cmd1)
    except:
        continue

    bpz  = final_folder + os.path.basename(catalog)[:-3]
    bpz += 'Mabs%s.bpz'%(filter_nick_name)
    if not os.path.exists(bpz):
       cmd  = 'python %sbpz.py %s '%(root2bpz,catalog)
       cmd += '-COLUMNS %s -OUTPUT %s -CHECK yes '%(columns,bpz)
       cmd += '-SIGMA_EXPECTED 0.015 -INTERP 15 -FLUX_COMPARISON yes '
       cmd += '-SPECTRA COSMOSeB11new_recal.list -ZMIN 0.005 -PRIOR SM '
       if filter_nick_name == 'rSDSS':
          cmd += '-DZ 0.001 -ZMAX 1.0 -HDF5 no -STELLAR_MASS yes -ONLY_TYPE no '
       else:
          cmd += '-DZ 0.001 -ZMAX 1.0 -HDF5 no -STELLAR_MASS no -ONLY_TYPE no '
       cmd += '-ABSOLUTE_MAGNITUDE yes ABSOLUTE_MAGNITUDE_FILTER %s '%(filters[ii][:-4])
       cmd += '-ZP_ERRORS "0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02,0.02" '
       print cmd
       os.system(cmd)



filtros = ['uJAVA','F0378W','F0395W','F0410W','F0430W','gSDSS',
          'F0515W','rSDSS','F0660W','iSDSS','F0861W','zSDSS']


ra,dec = U.get_data(catalog,(2,3))
MA1 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[0],filtros[0]),6)
MA2 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[1],filtros[1]),6)
MA3 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[2],filtros[2]),6)
MA4 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[3],filtros[3]),6)
MA5 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[4],filtros[4]),6)
MA6 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[5],filtros[5]),6)
MA7 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[6],filtros[6]),6)
MA8 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[7],filtros[7]),7)
MA9 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[8],filtros[8]),6)
MA10 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[9],filtros[9]),6)
MA11 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[10],filtros[10]),6)
MA12 = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[11],filtros[11]),6)
Mst,zs,Tb,Odds,mo = U.get_data(root2+'%s/Master.Mabs%s.bpz'%(filtros[7],filtros[7]),(6,11,4,5,12))

heada = '# RA Dec M_uJAVA M_J0378 M_J0395 M_J0410 M_J0430 M_gSDSS M_J0515 M_rSDSS '
heada += 'M_J0660 M_iSDSS M_J0861 M_zSDSS StellMass Specz SpT Odds rSDSS '
new_file = root2+'SPLUS_Master_MagAbs.cat'
U.put_data(new_file,(ra,dec,MA1,MA2,MA3,MA4,MA5,MA6,MA7,MA8,MA9,MA10,MA11,MA12,Mst,zs,Tb,Odds,mo),heada)
