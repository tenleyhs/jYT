'''
plots a 1D profile from FLASH and compares to the MESA input profile sm.dat

run like:
python jYT/initial_profile.py --cluster=fend --run=m1.0_p16_b3.0
'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yt
import os
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='hyades')
parser.add_argument('--run', help='run to plot, e.g., m1.0_p16_b3.0', type=str)
args = parser.parse_args()

if args.cluster == 'hyades':
	execfile('/pfs/lawsmith/jYT/my_settings.py')
	clusterdir = '/pfs/lawsmith/FLASH4.3/runs/'
	savepath = '/pfs/lawsmith/results/initial_profile/'
elif args.cluster == 'fend':
	execfile('/groups/dark/lawsmith/jYT/my_settings.py')
	clusterdir = '/groups/dark/lawsmith/FLASH4.3_copy/runs/'
	savepath = '/groups/dark/lawsmith/results/initial_profile/'
elif args.cluster == 'pfe':
	clusterdir = '/nobackup/jlawsmit/'
	savepath = '/nobackup/jlawsmit/results/initial_profile/'

plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

f = np.loadtxt(clusterdir + args.run + '/sm.dat', skiprows=1)
R = f[:,0]
rho = f[:,1]
T = f[:,2]
P = f[:,3]

chks = ['0000', '0005']

for chk in chks:
	ds = yt.load(clusterdir + args.run + '/multitidal_hdf5_chk_' + chk)
	sp = ds.sphere("c", (1.1*max(R), "cm"))

	plot = yt.ProfilePlot(sp, "radius", ["density"],
		weight_field="cell_mass", n_bins=2000, accumulation=False, x_log=False, y_log={'density':False})

	profile = plot.profiles[0]
	R_FLASH =  profile.x
	rho_FLASH = profile["density"]

	# account for zeros in FLASH array due to binning
	R_FLASH = R_FLASH[np.where(rho_FLASH > 0.)[0]]
	rho_FLASH = rho_FLASH[np.where(rho_FLASH > 0.)[0]]

	fig, ax = plt.subplots()
	ax.plot(R_FLASH/R_sun, np.log10(rho_FLASH), color='C0', label='FLASH')
	ax.plot(R/R_sun, np.log10(rho), color='C1', label='MESA')
	ax.set_xlabel(r'$R/R_\odot$')
	ax.set_ylabel(r'$\log\ \rho$')

	ax.legend()
	fig.tight_layout(pad=0.3)
	plt.savefig(savepath + args.run.replace(".","_") + '_' + chk + '.pdf')
