__author__ = 'albertomolino'

import os,sys
sys.path.append('/Users/albertomolino/doctorado/photo/programas/')
import numpy as N
import useful as U

root_to_files = '/Users/albertomolino/Desktop/'
fields_names_file = root_to_files+'stripe82_reduced_fields_12bands.cat'



# It creates a new file where to save the info.
new_filename = root_to_files+'new_tiles.list'
new_file = open(new_filename,'w')