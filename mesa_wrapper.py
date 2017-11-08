import argparse

parser = argparse.ArgumentParser(description='generate a file to input directly into FLASH from a MESA profile')
parser.add_argument('--indata', '-i', dest='indata', help='path and filename of MESA profile we want to prune', type=str)
parser.add_argument('--outdata', '-o', dest='outdata', help='path and filename of the output file we will feed to FLASH', type=str, default='sm.dat')
parser.add_argument('--elements', '-e', dest='els', help='elements we wish to track', type=str, default=['h1', 'he4', 'c12', 'o16'], nargs='+')
parser.add_argument('--plot', '-p', dest='plot', help='whether to save a plot or not', action='store_true', default=False)
args = parser.parse_args()

if not (args.indata):
	parser.error('No indata supplied. Add "-i mesadata.dat", for example.')

import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
import numpy as np

R_sun = 6.96E10


prof = ascii.read(args.indata, header_start=4, data_start=5, format='basic', guess=False)
#print prof.colnames

# plot something if we want to
if args.plot:
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.plot(prof['logRho'],prof['logP'],'-', color='black', lw=2)
	fig.savefig('P_vs_Rho.png')

# mass fraction of each of the elements we wish to track
mf_els = np.empty([len(args.els), len(prof)])
for i in np.arange(0, len(args.els)):
	mf_els[i] = prof[args.els[i]]

# make sure the mass fractions sum to unity
diff_max = 0
for row in np.arange(0,len(prof)):
	mf_sum = 0
	for column in np.arange(0,len(args.els)):
		mf_sum += prof[args.els][row][column]
	diff = (1.0-mf_sum)
	if diff > diff_max:
		diff_max = diff
if diff_max > 0.0:
	pass
	#print "WARNING: mass fractions of the elements do not always sum to unity. Greatest deviation is %f" % diff_max


# convert from log and into [cgs]
arrays = np.append([R_sun*10**prof['logR'], 10**prof['logRho'], 10**prof['logT'], 10**prof['logP']], mf_els, axis=0)
# reverse of the order of each of the columns so that we're going from small R to big R
arrays = np.fliplr(arrays)
names = np.append(['R', 'Rho', 'T', 'P'], args.els, axis=0)
data = Table(list(arrays), names=names)
# don't write col headers
ascii.write(data, args.outdata, format='no_header')

# want two integers in the first line listing the number of rows followed by
# the number of columns in the table. Then, output the table just below that.
new_first_line = str(len(prof['logR'])) + ' ' + str(len(names)) + '\n'
with open(args.outdata, 'r') as original: data = original.read()
with open(args.outdata, 'w') as modified: modified.write(new_first_line + data)
