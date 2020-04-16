'''
plots a 1D profile from FLASH and compares to the MESA input profile sm.dat

run like:
python jYT/initial_profile.py --cluster=fend --run m1.0_p16_b3.0 ... --chk 5 
'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import yt
import os
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., fend/lux/pfe', type=str, default='fend')
parser.add_argument('--run', help='run to plot, e.g., m1.0_p16_b3.0', type=str, nargs='+')
parser.add_argument('--chk', help='checkpoints to plot, e.g., 0 5', type=str, nargs='+')
parser.add_argument('--linear', help='linear instead of log', action='store_true', default=False)
args = parser.parse_args()

if args.cluster == 'fend':
	#execfile('/groups/dark/lawsmith/jYT/my_settings.py')
	exec(open('/groups/dark/lawsmith/jYT/my_settings.py').read())
	clusterdir = '/storage/dark/lawsmith/FLASH4.3/runs/'
	savepath = '/groups/dark/lawsmith/results/initial_profile/'
elif args.cluster == 'lux':
	#execfile('/home/lawsmith/jYT/my_settings.py')
	exec(open('/home/lawsmith/jYT/my_settings.py').read())
	clusterdir = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/'
	savepath = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/results/initial_profile/'
elif args.cluster == 'pfe':
	clusterdir = '/nobackup/jlawsmit/'
	savepath = '/nobackup/jlawsmit/results/initial_profile/'

for run in args.run:

	f = np.loadtxt(clusterdir + run + '/sm.dat', skiprows=1)
	R = f[:,0]
	rho = f[:,1]
	T = f[:,2]
	P = f[:,3]

	chks = [x.zfill(4) for x in args.chk]
	for chk in chks:
		ds = yt.load(clusterdir + run + '/multitidal_hdf5_chk_' + chk)
		sp = ds.sphere("c", (1.1*max(R), "cm"))

		plot = yt.ProfilePlot(sp, "radius", ["density"],
			weight_field="cell_mass", n_bins=2000, accumulation=False, x_log=False, y_log={'density':False})

		profile = plot.profiles[0]
		R_FLASH =  profile.x
		rho_FLASH = profile["density"]

		# account for zeros in FLASH array due to binning
		R_FLASH = R_FLASH[np.where(rho_FLASH > 0.)[0]]
		rho_FLASH = rho_FLASH[np.where(rho_FLASH > 0.)[0]]

		# prepend first density to zero. this is the density of the first cell
		R_FLASH = np.insert(np.array(R_FLASH), 0, 0)
		rho_FLASH = np.insert(np.array(rho_FLASH), 0, rho_FLASH[0])

		fig, ax = plt.subplots()
		if not args.linear:
			ax.plot(R/R_sun, np.log10(rho), color='k', lw=3, label='MESA')
			ax.plot(R_FLASH/R_sun, np.log10(rho_FLASH), color='r', lw=3, label='FLASH')
			ax.set_ylabel(r'$\log\ \rho\ {\rm [g/cm^3]}$')
		else:
			ax.plot(R/R_sun, rho, color='k', lw=3, label='MESA')
			ax.plot(R_FLASH/R_sun, rho_FLASH, color='r', lw=3, label='FLASH')
			ax.set_ylabel(r'$\rho\ {\rm [g/cm^3]}$')

		ax.set_xlabel(r'$R/R_\odot$')
		ax.grid(which='both')
		ax.set_ylim(bottom=min(np.log10(rho))-0.5)
		#ax.set_title(run + '; chk_' + chk + '; ' + str(round(ds.current_time.v))+'s', fontsize=14)
		ax.text(0.6, 0.9, '_'.join(run.split('_')[:2]), transform=ax.transAxes)
		ax.legend(loc='lower left')
		fig.tight_layout(pad=0.3)
		if not args.linear:
			plt.savefig(savepath + run.replace(".","_") + '_' + chk + '.pdf')
		else:
			plt.savefig(savepath + run.replace(".","_") + '_' + chk + '_linear.pdf')

		plt.close()
