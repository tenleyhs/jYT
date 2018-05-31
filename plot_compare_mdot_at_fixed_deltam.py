import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import ascii
execfile('/pfs/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18


"""
output from plot_deltaM_vs_beta
p0 -0.5 1.248
p0 -0.25 1.464
p0 -0.1 1.713
p10 -0.5 1.521
p10 -0.25 1.889
p10 -0.1 2.298
p16 -0.5 1.785
p16 -0.25 2.36
p16 -0.1 3.004
"""

# TODO should make this automatically from file
betas = {}
betas = {'p1': {'-0.5': 1.248,
                '-0.25': 1.464,
                '-0.1': 1.713},
        'p10': {'-0.5': 1.521,
                '-0.25': 1.889,
                '-0.1': 2.298},
        'p16': {'-0.5': 1.785,
                '-0.25': 2.363,
                '-0.1': 3.004},
        }

ps = ['p1', 'p10', 'p16']
dms = [-0.5,-0.25,-0.1]

for dm in dms:
    fig, ax = plt.subplots()
    ax.text(0, 1, str(dm), transform=ax.transAxes)
    #text(0.5, 0.5, 'matplotlib', horizontalalignment='center',
    #      verticalalignment='center', transform=ax.transAxes)
    for p in ps:
        beta = betas[p][str(dm)]
        #print beta
        dmdtdir = '/pfs/lawsmith/FLASH4.3/runs/results/dmdts/interpolated_m1_0_'+p+'_allbeta_0040_ev/'
        beta_files = os.listdir(dmdtdir)
        beta_files_floats = [float(b[:-4]) for b in beta_files]
        idx = (np.abs(np.array(beta_files_floats) - beta)).argmin()

        t, dmdt = np.loadtxt(dmdtdir + beta_files[idx], unpack=True)
        ax.plot(t, dmdt, label=p)

    #ax.set_xlim(-2, 0.25)
    #ax.set_ylim(-2, 1)
    ax.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
    ax.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
    ax.legend(loc=1)
    fig.tight_layout()
    fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/compare_mdot_at_fixed_deltam_'
                + str(dm).replace(".","_") + '.pdf')
    plt.close('all')
