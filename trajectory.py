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


fig,ax=plt.subplots()

#ax.plot(f['[02]posx'][::2], f['[03]posy'][::2], rasterized=True)
#ax.plot(f['[02]posx'][1::2], f['[03]posy'][1::2], rasterized=True)

# TODO for some reason this doesn't work with ::2, somethign to do with first column.
# A: aha! has to do with restarts probably
posx1 = np.array([float(x) for x in f['[02]posx'][::2]])/R_SUN_CGS
posx2 = np.array([float(x) for x in f['[02]posx'][1::2]])/R_SUN_CGS

posy1 = np.array([float(x) for x in f['[03]posy'][::2]])/R_SUN_CGS
posy2 = np.array([float(x) for x in f['[03]posy'][1::2]])/R_SUN_CGS

xcenter = np.mean(posx1)
ycenter = np.mean(posy1)

ax.plot(posx1 - xcenter, posy1 - ycenter, c='r', lw=3)
ax.plot(posx2 - xcenter, posy2 - ycenter, c='k', lw=3)

ax.set_xlabel(r'$x\ {\rm [R_\odot]}$')
ax.set_ylabel(r'$y\ {\rm [R_\odot]}$')
ax.grid()
plt.savefig(savepath +args.run+'/trajectory.pdf', bbox_inches='tight')

