"""
make all the histograms with a single python command inside the directory
there are nicer and maybe faster ways to paralellize this from googling
"""

import os

chk = '0090'
bins = '10000'

os.system('python /pfs/lawsmith/jYT/histogram.py multitidal_hdf5_chk_'+chk+' -v bhbound -n '+bins+' -ev selfbound -et 0. -ey min --silent -fp b'+bins+'_ev >> outputfile &')
os.system('python /pfs/lawsmith/jYT/histogram_h1.py multitidal_hdf5_chk_'+chk+' -v bhbound -n '+bins+' -ev selfbound -et 0. -ey min --silent -fp b'+bins+'_h1 >> outputfile &')
os.system('python /pfs/lawsmith/jYT/histogram_he4.py multitidal_hdf5_chk_'+chk+' -v bhbound -n '+bins+' -ev selfbound -et 0. -ey min --silent -fp b'+bins+'_he4 >> outputfile &')
os.system('python /pfs/lawsmith/jYT/histogram_o16.py multitidal_hdf5_chk_'+chk+' -v bhbound -n '+bins+' -ev selfbound -et 0. -ey min --silent -fp b'+bins+'_o16 >> outputfile &')
os.system('python /pfs/lawsmith/jYT/histogram_c12.py multitidal_hdf5_chk_'+chk+' -v bhbound -n '+bins+' -ev selfbound -et 0. -ey min --silent -fp b'+bins+'_c12 >> outputfile &')
os.system('python /pfs/lawsmith/jYT/histogram_ne20.py multitidal_hdf5_chk_'+chk+' -v bhbound -n '+bins+' -ev selfbound -et 0. -ey min --silent -fp b'+bins+'_ne20 >> outputfile &')
os.system('python /pfs/lawsmith/jYT/histogram_n14.py multitidal_hdf5_chk_'+chk+' -v bhbound -n '+bins+' -ev selfbound -et 0. -ey min --silent -fp b'+bins+'_n14 >> outputfile &')
