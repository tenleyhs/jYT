"""
plot Mass in annulus inside secondary orbit, but outside core radius

don't run in parallel, just python ~/jYT/M_inside_secondary_outside_core.py ...
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
import yt
import os
import argparse
import astropy.constants as c
R_SUN_CGS = c.R_sun.cgs.value
M_SUN_CGS = c.M_sun.cgs.value
G = c.G.cgs.value

### TODO hardcoded in for G18p201 right now ###

R_CORE = 0.8
M_CORE = 4.359649747506503
M_TOTAL = 4.56318599
M_ENV = M_TOTAL-M_CORE
XCENTER = 140.e10
YCENTER = 140.e10
ZCENTER = 140.e10
sink_softening_radius = 1.16e10

### ###


parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='fend')
parser.add_argument('--run', help='run to plot, e.g., m1.0_p16_b3.0', type=str)
parser.add_argument('--chkplt', help='plot from chks or plts', type=str, default='chk')
parser.add_argument('--files', help='which files to plot, e.g. [0-9][0-9][0-9]0 or *', type=str, default='*')
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


def _unbound_mass(field, data):
    e = 0.5*(data['velocity_magnitude'].d)**2 + data['gpot'].d
    
    factor = np.sign(e)

    return factor*data['cell_mass'].d/M_SUN_CGS

yt.add_field(("gas","unbound_mass"), display_name=r'$M_{\rm ej}$', function=_unbound_mass, take_log=False, units=None, force_override=True)


if args.chkplt == 'chk':
    LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_chk_' + args.files
elif args.chkplt == 'plt':
    LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_plt_cnt_' + args.files

yt.enable_parallelism()
ts = yt.DatasetSeries(LOAD_FILES)

t_array = []
m_array = []

f = ascii.read(clusterdir + args.run + '/' + args.input)

time_array = f['[01]time'][1::2]
radius_array = np.sqrt((f['[02]posx'][1::2]-XCENTER)**2 +(f['[03]posy'][1::2]-YCENTER)**2 + (f['[04]posz'][1::2]-ZCENTER)**2)

for ds in ts.piter():

    t = float(ds.current_time.d)
    print(t)
    t_array.append(t)

    index = (np.abs(time_array - t)).argmin()
    r = radius_array[index]

    ad = ds.all_data()

    # including softening radius
    #string1 = "obj['radius'] < " + str(r-sink_softening_radius)
    string1 = "obj['radius'] < 20e10"
    inside_ad = ad.cut_region([string1])
    #string2 = "obj['radius'] > " + str(R_CORE*R_SUN_CGS)
    string2 = "obj['radius'] > 5.6e10"
    outside_and_inside_ad = inside_ad.cut_region([string2])
    
    #inside_ad = ad.include_below('radius', r-sink_softening_radius)
    #outside_and_inside_ad = inside_ad.include_above('radius', R_CORE*R_SUN_CGS)

    mass = outside_and_inside_ad.quantities.total_quantity(["cell_mass"]).d/M_SUN_CGS
    print(mass)
    m_array.append(mass)



fig, ax = plt.subplots()

ax.plot(np.array(t_array)/3600., np.array(m_array), lw=5, c='k', label=r'$M\ (R_{\rm core} < r < R_{\rm secondary})$')
#ax.axhline(M_ENV, c='k', ls='--', lw=5, label=r'$M_{\rm envelope}$')
#ax.axhline(M_CORE, c='k', ls=':', lw=6, label=r'$M_{\rm core}$')

ax.set_xlim(left=0)
ax.set_xlabel(r'$t\ {\rm [h]}$')
ax.set_ylabel(r'$M\ {\rm [M_\odot]}$')
ax.legend()
ax.grid()
fig.tight_layout(pad=0.3)
plt.savefig(savepath +args.run+'/M_inside_secondary_outside_core.pdf')