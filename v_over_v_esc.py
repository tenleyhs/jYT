"""
plot velocity magnitude over escape velocity
"""

# start with escape velocity of core. did this in jupyter notebook
#for G17p55, core is 0.31 Rsun, 3.44 Msun. vesc=205728746 cm/s
#for G18p201, core is 0.8 Rsun, 4.36 Msun. vesc=144240990 cm/s

# todo --chkplt=plt doesn't work yet bc not saving gpot to plt files. could change this.

# todo add x at position of secondary. just need to read in sinks_evol.dat and take nearest time

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


def _v_over_v_esc(field, data):
    v_esc = np.sqrt(2.*np.abs(data['gpot'].d))
    return data['velocity_magnitude'].d/v_esc
    

yt.add_field(("gas","v_over_v_esc"), display_name='|v|/v_{esc,local}', function=_v_over_v_esc, take_log=False, units=None, force_override=True)

if args.chkplt == 'chk':
    LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_chk_' + args.files
elif args.chkplt == 'plt':
    LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_plt_cnt_' + args.files

yt.enable_parallelism()
ts = yt.DatasetSeries(LOAD_FILES)

for ds in ts.piter():

    #ds = yt.load(clusterdir + args.run + '/' ......

    #ad = ds.all_data()
    #cut_ad = ad.cut_region(["obj['v_over_v_esc'] >= 1."])
    #s = yt.SlicePlot(ds, 'z', 'v_over_v_esc', data_source=cut_ad)

    s = yt.SlicePlot(ds, 'z', 'v_over_v_esc')

    # doesn't work
    #s.set_background_color('v_over_v_esc', 'white')

    s.annotate_contour("v_over_v_esc", ncont=1, clim=(1,1), label=False,
                        plot_args={"colors":"white","linewidths":2,"linestyles":"-"},
                        text_args={"size":48})

    s.annotate_timestamp()
    s.annotate_scale(unit='rsun')

    if args.vel:
        s.annotate_velocity()

    if args.line_integral:
        s.annotate_line_integral_convolution('velx', 'vely')#, lim=(0.5,0.65))

    if argsr.streamlines:
        s.annotate_streamlines('velx', 'vely')

    if args.cell_edges:
        s.annotate_cell_edges()
    
    if args.grids:
        s.annotate_grids()

    if args.vel:
        s.annotate_velocity()
        s.save(savepath + args.run + '/v_over_v_esc_vel/')

    elif args.streamlines:
        s.annotate_streamlines('velx', 'vely')
        s.save(savepath + args.run + '/v_over_v_esc_streamlines/')

    elif args.line_integral:
        s.annotate_line_integral_convolution('velx', 'vely')#, lim=(0.5,0.65))
        s.save(savepath + args.run + '/v_over_v_esc_line_integral/')

    else:
        s.save(savepath + args.run + '/v_over_v_esc/')