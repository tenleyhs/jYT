"""
plot E_bind (= ekin + gpot) vs time
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
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
E_BIND_CORE = G*M_CORE*M_SUN_CGS/(R_CORE*R_SUN_CGS)

### ###


parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='fend')
parser.add_argument('--run', help='run to plot, e.g., m1.0_p16_b3.0', type=str)
parser.add_argument('--chkplt', help='plot from chks or plts', type=str, default='chk')
parser.add_argument('--files', help='which files to plot, e.g. [0-9][0-9][0-9]0 or *', type=str, default='*')
parser.add_argument('--vel', help='annotate_velocity?', type=bool, default=False)
parser.add_argument('--streamlines', help='annotate_streamlines?', type=bool, default=False)
parser.add_argument('--line_integral', help='annotate_line_integral_convolution?', type=bool, default=False)
parser.add_argument('--cell_edges', help='overplot cell edges?', type=bool, default=False)
parser.add_argument('--grids', help='overplot grid?', type=bool, default=False)


args = parser.parse_args()

if args.cluster == 'fend':
    exec(open('/groups/dark/lawsmith/jYT/my_settings.py').read())
    clusterdir = '/storage/dark/lawsmith/FLASH4.3/runs/'
    savepath = '/groups/dark/lawsmith/fend_slices/'
elif args.cluster == 'lux':
    exec(open('/home/lawsmith/jYT/my_settings.py').read())
    clusterdir = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/'
    savepath = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/lux_slices/'


def _ekin_plus_gpot(field, data):
    return 0.5*(data['velocity_magnitude'].d)**2 + data['gpot'].d

yt.add_field(("gas","ekin_plus_gpot"), display_name='v^2/2 + gpot', function=_ekin_plus_gpot, take_log=False, units=None, force_override=True)


if args.chkplt == 'chk':
    LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_chk_' + args.files
elif args.chkplt == 'plt':
    LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_plt_cnt_' + args.files

yt.enable_parallelism()
ts = yt.DatasetSeries(LOAD_FILES)

t_array = []
e_array = []

for ds in ts.piter():

    t_array.append(float(ds.current_time.d))

    ad = ds.all_data()

    outside_core = ad.cut_region(['obj["radius"] > 5.564e10']) #0.8 rsun

    e = outside_core.quantities.total_quantity(["ekin_plus_gpot"]).d

    e_array.append(e)


fig, ax = plt.subplots()

ax.plot(np.array(t_array)/3600., np.array(e_array), lw=5, label=r'$e_{\rm bind} = v^2/2 + GM_{\rm enc}/r$', c='k')
#ax.axhline(E_BIND_CORE, ls='--', lw=5, label=r'$GM_{\rm core}/R_{\rm core}$')

ax.set_xlabel(r'$t\ {\rm [h]}$')
ax.set_ylabel(r'$e\ {\rm [erg/g]}$')
ax.legend()
ax.grid()
fig.tight_layout(pad=0.3)
plt.savefig(savepath +args.run+'/E_bind_vs_time.pdf')