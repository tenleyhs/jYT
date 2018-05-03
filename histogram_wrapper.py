"""
make all the histograms with a single python command inside the directory
there are nicer and maybe faster ways to paralellize this from googling
"""

import os

os.system('python /pfs/lawsmith/jYT/histogram.py multitidal_hdf5_chk_0040 -v bhbound -n 10000 -ev selfbound -et 0. -ey min --silent -fp b10000_ev >> outputfile_ev &')
os.system('python /pfs/lawsmith/jYT/histogram_h1.py multitidal_hdf5_chk_0040 -v bhbound -n 10000 -ev selfbound -et 0. -ey min --silent -fp b10000_h1 >> outputfile_h1 &')
os.system('python /pfs/lawsmith/jYT/histogram_he4.py multitidal_hdf5_chk_0040 -v bhbound -n 10000 -ev selfbound -et 0. -ey min --silent -fp b10000_he4 >> outputfile_he4 &')
os.system('python /pfs/lawsmith/jYT/histogram_o16.py multitidal_hdf5_chk_0040 -v bhbound -n 10000 -ev selfbound -et 0. -ey min --silent -fp b10000_o16 >> outputfile_o16 &')
os.system('python /pfs/lawsmith/jYT/histogram_c12.py multitidal_hdf5_chk_0040 -v bhbound -n 10000 -ev selfbound -et 0. -ey min --silent -fp b10000_c12 >> outputfile_c12 &')
os.system('python /pfs/lawsmith/jYT/histogram_ne20.py multitidal_hdf5_chk_0040 -v bhbound -n 10000 -ev selfbound -et 0. -ey min --silent -fp b10000_ne20 >> outputfile_ne20 &')
os.system('python /pfs/lawsmith/jYT/histogram_n14.py multitidal_hdf5_chk_0040 -v bhbound -n 10000 -ev selfbound -et 0. -ey min --silent -fp b10000_n14 >> outputfile_n14 &')
