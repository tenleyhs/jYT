"""
plot separation a of two sinks vs time
"""

import matplotlib
matplotlib.use('Agg')
import argparse
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import astropy.constants as c
R_SUN_CGS = c.R_sun.cgs.value

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='fend')
parser.add_argument('--run', help='run to plot, e.g., m1.0_p16_b3.0', type=str)
parser.add_argument('--input', help='option when doesnt work for restarts in sinks_evol.dat', type=str, default='sinks_evol.dat')

args = parser.parse_args()

if args.cluster == 'fend':
    exec(open('/groups/dark/lawsmith/jYT/my_settings.py').read())
    clusterdir = '/storage/dark/lawsmith/FLASH4.3/runs/'
    savepath = '/groups/dark/lawsmith/fend_slices/'
elif args.cluster == 'lux':
    exec(open('/home/lawsmith/jYT/my_settings.py').read())
    clusterdir = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/'
    savepath = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/lux_slices/'


f = ascii.read(clusterdir + args.run + '/' + args.input)

dist = np.sqrt((f['[02]posx'][::2] - f['[02]posx'][1::2])**2 \
               + (f['[03]posy'][::2] - f['[03]posy'][1::2])**2 \
               + (f['[04]posz'][::2] - f['[04]posz'][1::2])**2)


fig,ax=plt.subplots()

ax.plot(f['[01]time'][::2]/3600., dist/R_SUN_CGS, lw=3, c='k')

ax.set_ylim(bottom=0)
ax.set_xlabel(r'$t\ [h]$')
ax.set_ylabel(r'$a\ [R_\odot]$')
ax.grid()
fig.tight_layout(pad=0.3)
plt.savefig(savepath +args.run+'/a_vs_t.pdf')

ax.set_yscale('log')
ax.set_ylim(1e-1, 10)
ax.grid()
fig.tight_layout(pad=0.3)
plt.savefig(savepath +args.run+'/a_vs_t_log.pdf')

