"""
plot M_ej vs time

don't run in parallel, just python ~/jYT/M_ej_vs_time.py ...
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

M_CORE = 4.359649747506503
M_TOTAL = 4.56318599
M_ENV = M_TOTAL-M_CORE

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
unbound_m_array = []
bound_m_array = []

for ds in ts.piter():

    t_array.append(float(ds.current_time.d))

    #s = yt.SlicePlot(ds, 'z', 'unbound_mass')
    #s.save(savepath + args.run + '/unbound_mass/')

    ad = ds.all_data()

    pos_ad = ad.cut_region(['obj["unbound_mass"] >= 0.'])
    neg_ad = ad.cut_region(['obj["unbound_mass"] < 0.'])
    unbound_mass = pos_ad.quantities.total_quantity(["unbound_mass"]).d
    bound_mass = neg_ad.quantities.total_quantity(["unbound_mass"]).d
    unbound_m_array.append(float(unbound_mass))
    bound_m_array.append(float(bound_mass))




fig, ax = plt.subplots()

ax.plot(np.array(t_array)/3600., np.array(unbound_m_array), lw=5, label=r'$M_{\rm ej}$')
ax.axhline(M_ENV, c='k', ls='--', lw=5, label=r'$M_{\rm envelope}$')
ax.plot(np.array(t_array)/3600., -np.array(bound_m_array), lw=5, label=r'$M_{\rm bound}$')
ax.axhline(M_CORE, c='k', ls=':', lw=6, label=r'$M_{\rm core}$')

ax.set_xlim(left=0)
ax.set_xlabel(r'$t\ {\rm [h]}$')
ax.set_ylabel(r'$M\ {\rm [M_\odot]}$')
ax.legend()
ax.grid()
fig.tight_layout(pad=0.3)
plt.savefig(savepath +args.run+'/M_ej_vs_time.pdf')

fig, ax = plt.subplots()
ax.plot(np.array(t_array)/3600., np.array(unbound_m_array)/M_ENV, lw=5, label=r'$M_{\rm ej}/M_{\rm envelope}$')
ax.plot(np.array(t_array)/3600., -np.array(bound_m_array)/M_CORE, lw=5, label=r'$M_{\rm bound}/M_{\rm core}$')

ax.set_xlim(left=0)
ax.set_xlabel(r'$t\ {\rm [h]}$')
ax.set_ylabel('fraction')
ax.legend()
ax.grid()
fig.tight_layout(pad=0.3)
plt.savefig(savepath +args.run+'/M_ej_vs_time_fractional.pdf')