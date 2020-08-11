"""
make all the histograms with a single python command inside the directory
there are nicer and maybe faster ways to paralellize this from googling

run inside diretory like:
python ~/jYT/histogram_wrapper.py --cluster=fend --chk=0080 --bins=10000
"""
import os
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., fend/pfe', type=str, default='fend')
parser.add_argument('--chk', help='checkpoint to plot, e.g., 0080', type=str)
parser.add_argument('--bins', help='number of bins in histogram, e.g., 10000', type=str, default='10000')
args = parser.parse_args()

if args.cluster == 'fend': 
    path = '/groups/dark/lawsmith/jYT/'
elif args.cluster == 'lux':
    path = '/home/lawsmith/jYT/'

# TODO should probably just submit these to cluster. Or have a switch to do on head nodes or cluster.

# not sure about > vs >>. Both > and >> output the progress to the terminal, not outputfile. But then print bounds of histogram in outputfile.
#os.system('python '+path+'histogram.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_ev >> outputfile &')

os.system('python '+path+'histogram_h1.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_h1 >> outputfile &')
os.system('python '+path+'histogram_he4.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_he4 >> outputfile &')
os.system('python '+path+'histogram_o16.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_o16 >> outputfile &')
os.system('python '+path+'histogram_c12.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_c12 >> outputfile &')
os.system('python '+path+'histogram_ne20.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_ne20 >> outputfile &')
os.system('python '+path+'histogram_n14.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_n14 >> outputfile &')

#os.system('python '+path+'histogram_he3.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_he3 >> outputfile &')
#os.system('python '+path+'histogram_c13.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_c13 >> outputfile &')
#os.system('python '+path+'histogram_na23.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_na23 >> outputfile &')
#os.system('python '+path+'histogram_mg24.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_mg24 >> outputfile &')
#os.system('python '+path+'histogram_al27.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_al27 >> outputfile &')
#os.system('python '+path+'histogram_si28.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_si28 >> outputfile &')
#os.system('python '+path+'histogram_s34.py multitidal_hdf5_chk_'+args.chk+' -v bhbound -n '+args.bins+' -ev selfbound -et 0. -ey min --silent -fp b'+args.bins+'_s34 >> outputfile &')
