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
p1 -0.5 1.2605
p1 -0.25 1.4769
p1 -0.1 1.6933
p10 -0.5 1.4929
p10 -0.25 1.8116
p10 -0.1 2.1603
p16 -0.5 1.7034
p16 -0.25 2.1603
p16 -0.1 2.7074
"""

# TODO should make this automatically from file
betas = {}
betas = {'p1': {'-0.5': 1.2605,
                '-0.25': 1.4769,
                '-0.1': 1.6933},
        'p10': {'-0.5': 1.4929,
                '-0.25': 1.8116,
                '-0.1': 2.1603},
        'p16': {'-0.5': 1.7034,
                '-0.25': 2.1603,
                '-0.1': 2.7074},
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
        dmdtdir = '/pfs/lawsmith/FLASH4.3/runs/results/dmdts/interpolated_m1_0_'+p+'_allbeta_lastchk_ev/'
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
    fig.tight_layout(pad=0.3)
    fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/paper/compare_mdot_at_fixed_deltam_'
                + str(dm).replace(".","_") + '.pdf')
    plt.close('all')
