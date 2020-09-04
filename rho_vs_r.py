"""
plot separation rho vs r,t

don't run in parallel, just:
python ~/jYT/rho_vs_r.py --run=mm_G18p201_15 --files=00*[0,5]
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

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='fend')
parser.add_argument('--run', help='run to plot, e.g., m1.0_p16_b3.0', type=str)
parser.add_argument('--input', help='option when doesnt work for restarts in sinks_evol.dat', type=str, default='sinks_evol.dat')
parser.add_argument('--chkplt', help='plot from chks or plts', type=str, default='chk')
parser.add_argument('--files', help='which files to plot, e.g. [0-9][0-9][0-9]0 or *', type=str, default='*')

args = parser.parse_args()

if args.cluster == 'fend':
    exec(open('/groups/dark/lawsmith/jYT/my_settings.py').read())
    clusterdir = '/storage/dark/lawsmith/FLASH4.3/runs/'
    savepath = '/groups/dark/lawsmith/fend_slices/'
elif args.cluster == 'lux':
    exec(open('/home/lawsmith/jYT/my_settings.py').read())
    clusterdir = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/'
    savepath = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/lux_slices/'


if args.chkplt == 'chk':
    LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_chk_' + args.files
elif args.chkplt == 'plt':
    LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_plt_cnt_' + args.files

yt.enable_parallelism()
ts = yt.DatasetSeries(LOAD_FILES)

fig, ax = plt.subplots()

for ds in ts.piter():

    sp0 = ds.sphere(ds.domain_center, (1, "Mpc"))
    rp0 = yt.create_profile(sp0, 'radius', 'dens', logs = {'radius': False})
    ax.plot(rp0.x.value/R_SUN_CGS, rp0["dens"].value, label=str(round(float(ds.current_time.d)/3600.,1))+'hr', lw=5)
    #ad = ds.all_data()
    #profiles.append(yt.create_profile(ad, ["radius"], fields=["dens"], weight_field=None))
    #labels.append("t = %.2f" % ds.current_time)


ax.set_yscale('log')
ax.set_xlim(0, 15)
ax.set_ylim(bottom=1e-8)
ax.grid()
ax.legend()
ax.set_xlabel(r'$R\ {\rm [R_\odot]}$')
ax.set_ylabel(r'$\rho\ {\rm [g/cm^3]}$')
fig.savefig(savepath + args.run + '/rho_vs_r.pdf', bbox_inches='tight')

ax.set_xscale('log')
ax.set_xlim(0.2, 20)
ax.set_xlabel(r'$R\ {\rm [R_\odot]}$')
fig.savefig(savepath + args.run + '/rho_vs_r_log.pdf', bbox_inches='tight')

#p = yt.ProfilePlot.from_profiles(profiles, labels=labels)
#p.set_line_property("linewidth", "3")
#p.set_log("radius", False)
#p.save(savepath + args.run + '/')