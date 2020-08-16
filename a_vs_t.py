"""
plot separation a of two sinks vs time
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import astropy.constants as c
R_SUN_CGS = c.R_sun.cgs.value

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='fend')
parser.add_argument('--run', help='run to plot, e.g., m1.0_p16_b3.0', type=str)

args = parser.parse_args()

if args.cluster == 'fend':
    exec(open('/groups/dark/lawsmith/jYT/my_settings.py').read())
    clusterdir = '/storage/dark/lawsmith/FLASH4.3/runs/'
    savepath = '/groups/dark/lawsmith/fend_slices/'
elif args.cluster == 'lux':
    exec(open('/home/lawsmith/jYT/my_settings.py').read())
    clusterdir = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/'
    savepath = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/lux_slices/'



f = ascii.read(clusterdir + args.run + '/sinks_evol.dat')

dist = np.sqrt((f['[02]posx'][::2] - f['[02]posx'][1::2])**2 \
               + (f['[03]posy'][::2] - f['[03]posy'][1::2])**2 \
               + (f['[04]posz'][::2] - f['[04]posz'][1::2])**2)


fig,ax=plt.subplots()

ax.plot(f['[01]time'][::2], dist/R_SUN_CGS)

ax.set_xlabel(r'$t~[s]$')
ax.set_ylabel(r'$a~[R_\odot]$')
ax.grid(which='both')
fig.tight_layout(pad=0.3)
plt.savefig(savepath +args.run+'/a_vs_t.pdf')

