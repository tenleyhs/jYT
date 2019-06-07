"""
make all the histograms with a single python command inside the directory
there are nicer and maybe faster ways to paralellize this from googling

run inside diretory like:
python ../../../jYT/histogram_wrapper.py --cluster=fend --chk=0080 --bins=10000
"""
import os
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='fend')
parser.add_argument('--chk', help='checkpoint to plot, e.g., 0080', type=str)
parser.add_argument('--bins', help='number of bins in histogram, e.g., 10000', type=str, default='10000')
args = parser.parse_args()

if args.cluster == 'hyades': path = '/pfs/lawsmith/jYT/'
elif args.cluster == 'fend': path = '/groups/dark/lawsmith/jYT/'

#os.system('python '+path+'histogram.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_ev >> outputfile &')

os.system('python '+path+'histogram_h1.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_h1 >> outputfile &')
os.system('python '+path+'histogram_he4.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_he4 >> outputfile &')
os.system('python '+path+'histogram_o16.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_o16 >> outputfile &')
os.system('python '+path+'histogram_c12.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_c12 >> outputfile &')
os.system('python '+path+'histogram_ne20.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_ne20 >> outputfile &')
os.system('python '+path+'histogram_n14.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_n14 >> outputfile &')
