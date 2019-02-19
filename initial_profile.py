'''
plots a 1D profile from FLASH and compares to the MESA input profile sm.dat
'''

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt

execfile('my_settings.py')

PROF_PATH = '/pfs/lawsmith/FLASH4.3/runs_loadprofile/test/sm.dat'
CHK_PATH = '/pfs/lawsmith/FLASH4.3/runs/test/multitidal_hdf5_chk_0000'
SAVE_PNG_PATH = '/pfs/lawsmith/FLASH4.3/runs/test/'
lw = 1.5

f = np.loadtxt(PROF_PATH, skiprows=1)
R = f[:,0]
rho = f[:,1]
T = f[:,2]
P = f[:,3]

ds = yt.load(CHK_PATH)
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
ax.plot(R_FLASH/R_sun, np.log10(rho_FLASH), color=red, lw=lw, label='FLASH')
ax.plot(R/R_sun, np.log10(rho), color=blue, lw=lw, label='MESA')
ax.set_xlabel(r'$R/R_\odot$')
ax.set_ylabel(r'$\log\ \rho$')

ax.set_ylim(-8, -2)

ax.legend()
fig.tight_layout()
plt.savefig(SAVE_PNG_PATH + 'initial_profile.png')
