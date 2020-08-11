"""
butterfly plot Enrico described
"""

import yt
import colormaps as cmaps
import argparse, sys
import os

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., fend/lux', type=str, default='fend')
parser.add_argument('--t', help='tdyn after peri, e.g., 1/2/3', type=int, default=3)
args = parser.parse_args()

plots = [
#    ['m0.3_p1_b0.6_300k',   23, 0.2814, 10**2.0433069661766980, r'$0.3M_\odot,0Gyr,\beta=0.6$', 'f'],
#    ['m0.3_p1_b0.7_300k',   23, 0.2814, 10**2.0433069661766980, r'$0.3M_\odot,0Gyr,\beta=0.8$', 'f'],
#    ['m0.3_p1_b0.8_300k',   22, 0.2814, 10**2.0433069661766980, r'$0.3M_\odot,0Gyr,\beta=0.8$', 'f'],
#    ['m0.3_p1_b0.9_300k',   22, 0.2814, 10**2.0433069661766980, r'$0.3M_\odot,0Gyr,\beta=0.8$', 'f'],
#    ['m0.3_p1_b1.0_300k',   22, 0.2814, 10**2.0433069661766980, r'$0.3M_\odot,0Gyr,\beta=1.0 $', 'f'],
#    ['m0.3_p1_b2.0_300k',   21, 0.2814, 10**2.0433069661766980, r'$0.3M_\odot,0Gyr,\beta=2.0 $', 'l'],
##    
#    #['m0.5']
##
#    ['m0.7_p50_b0.8_300k',   22, 0.6793, 10**2.0600158274310392, r'$0.7M_\odot,10Gyr,\beta=1.0$',    'f'],
#    ['m0.7_p50_b1.0_300k',   22, 0.6793, 10**2.0600158274310392, r'$0.7M_\odot,10Gyr,\beta=1.0$',    'f'],
#    ['m0.7_p50_b1.15_300k',  22, 0.6793, 10**2.0600158274310392, r'$0.7M_\odot,10Gyr,\beta=1.0$',   'f'],
#    ['m0.7_p50_b1.3_300k',   22, 0.6793, 10**2.0600158274310392, r'$0.7M_\odot,10Gyr,\beta=1.0$',    'f'],
#    ['m0.7_p50_b1.5_300k',   21, 0.6793, 10**2.0600158274310392, r'$0.7M_\odot,10Gyr,\beta=1.0$',    'f'],
#    ['m0.7_p50_b1.7_300k',   21, 0.6793, 10**2.0600158274310392, r'$0.7M_\odot,10Gyr,\beta=1.0$',    'f'],
#    ['m0.7_p50_b3.6_300k',   21, 0.6793, 10**2.0600158274310392, r'$0.7M_\odot,10Gyr,\beta=1.0$',    'l'],
##
#    ['m1.0_p10_b1.0_300k',   22, 1.0455, 10**2.2301308164317524, r'$1M_\odot,4.8Gyr,\beta=1.0 $', 'f'],
#    ['m1.0_p10_b1.5_300k',   21, 1.0455, 10**2.2301308164317524, r'$1M_\odot,4.8Gyr,\beta=1.0 $', 'f'],
#    ['m1.0_p10_b2.0_300k',   21, 1.0455, 10**2.2301308164317524, r'$1M_\odot,4.8Gyr,\beta=1.0 $', 'f'],
#    ['m1.0_p10_b2.5_300k',   21, 1.0455, 10**2.2301308164317524, r'$1M_\odot,4.8Gyr,\beta=1.0 $', 'f'],
#    ['m1.0_p10_b3.0_300k',   20, 1.0455, 10**2.2301308164317524, r'$1M_\odot,4.8Gyr,\beta=1.0 $', 'l'],
#    ['m1.0_p10_b3.5_300k',   20, 1.0455, 10**2.2301308164317524, r'$1M_\odot,4.8Gyr,\beta=1.0 $', 'l'],
#    ['m1.0_p10_b4.9_300k',   20, 1.0455, 10**2.2301308164317524, r'$1M_\odot,4.8Gyr,\beta=1.0 $', 'l'],
##
#    ['m1.0_p16_b1.0_300k_lr16',  22, 1.2872, 10**2.6989544680278761, r'$1M_\odot,8.4Gyr,\beta=2.0$', 'f'],
#    ['m1.0_p16_b1.5_300k_lr16', 1071, 1.2872, 10**2.6989544680278761, r'$1M_\odot,8.4Gyr,\beta=2.0$', 'f'],
#    ['m1.0_p16_b2.0_300k_lr16',   21, 1.2872, 10**2.6989544680278761, r'$1M_\odot,8.4Gyr,\beta=2.0$', 'f'],  # run into problem with contours at 3 tdyn here
#    ['m1.0_p16_b3.0_300k_lr16',   20, 1.2872, 10**2.6989544680278761, r'$1M_\odot,8.4Gyr,\beta=3.0$', 'f'],
#    ['m1.0_p16_b4.0_300k_lr16',   20, 1.2872, 10**2.6989544680278761, r'$1M_\odot,8.4Gyr,\beta=3.0$', 'f'],
#    ['m1.0_p16_b5.0_300k_lr16',   20, 1.2872, 10**2.6989544680278761, r'$1M_\odot,8.4Gyr,\beta=5.0$', 'f'],
    ['m1.0_p16_b6.0_300k_lr16',   21, 1.2872, 10**2.6989544680278761, r'$1M_\odot,8.4Gyr,\beta=6.0$', 'f'],
##    
#    ['m3.0_p16_b1.5_300k_lr16',   21, 3.3192, 10**2.1421056939583099, r'$3M_\odot,2Gyr,\beta=2.0$', 'f'],
#    ['m3.0_p16_b2.0_300k_lr16',   21, 3.3192, 10**2.1421056939583099, r'$3M_\odot,2Gyr,\beta=2.0$', 'f'],
#    ['m3.0_p16_b3.0_300k_lr16',   21, 3.3192, 10**2.1421056939583099, r'$3M_\odot,2Gyr,\beta=3.0$', 'f'],
#    ['m3.0_p16_b4.0_300k_lr16',   22, 3.3192, 10**2.1421056939583099, r'$3M_\odot,2Gyr,\beta=3.0$', 'f'],
#    ['m3.0_p16_b5.0_300k_lr16',   20, 3.3192, 10**2.1421056939583099, r'$3M_\odot,2Gyr,\beta=5.0$', 'f'],
#    ['m3.0_p16_b7.0_300k_lr16',   20, 3.3192, 10**2.1421056939583099, r'$3M_\odot,2Gyr,\beta=7.0$', 'f'],
]

t = args.t
var = 'dens'
width = (10, 'rsun')
log_tf = True
set_zlim_tf = True
contour_tf = True
sz = 48

#LOAD_FILES = []

for plot in plots:
    if args.cluster == 'fend':
        if plot[-1] == 'l':
            continue
    if args.cluster == 'lux':
        if plot[-1] == 'f':
            continue

    if args.cluster == 'fend':
        LOAD_PATH = '/storage/dark/lawsmith/FLASH4.3/runs/'
    elif args.cluster == 'lux':
        LOAD_PATH = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/'
    
    # trying parallel
    #LOAD_FILES.append(LOAD_PATH+plot[0]+'/multitidal_hdf5_chk_'+str(plot[1]+t).zfill(4))
#yt.enable_parallelism()
#ts = yt.DatasetSeries(LOAD_FILES)
#for ds in ts.piter():

    ds = yt.load(LOAD_PATH+plot[0]+'/multitidal_hdf5_chk_'+str(plot[1]+t).zfill(4))

    width = (10*plot[2], 'rsun')
    s = yt.SlicePlot(ds, 'z', var, width=width)
    s.set_log(var, log_tf)
    if contour_tf:
        s.annotate_contour("density", ncont=3, clim=(1e-2,1e0), label=False,
                            plot_args={"colors":"white","linewidths":2,"linestyles":"-"},
                            text_args={"size":48})
    if set_zlim_tf:
        s.set_zlim(var, 1e-8*plot[3], plot[3])
    s.set_cmap(var, cmaps.viridis)
    #s.annotate_timestamp(time_unit='s')
    #if plot[0] == 'm1.0_p16_b2.0':
    #    s.annotate_scale(unit='rsun', text_args={"size":48})
    #s.set_axes_unit('rsun')
    #s.set_font({'size':24})
    #if plot[0] == 'm1.0_p1_b2.0':
    #    s.annotate_text((0.65, 0.92), plot[4], coord_system='axis',
    #        text_args={"size":60, "color":"white", "backgroundcolor":None})
    #else:
    #    s.annotate_text((0.6, 0.92), plot[4], coord_system='axis',
    #        text_args={"size":60, "color":"white", "backgroundcolor":None})
    s.hide_axes()
    s.hide_colorbar()
    s.set_figure_size(12)
    #s.set_buff_size(2000)
    if args.cluster == 'fend':
        SAVE_PATH = '/groups/dark/lawsmith/results/butterfly/1e-8/'+str(t)+'tdyn/'
    elif args.cluster == 'lux':
        SAVE_PATH = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/results/butterfly/1e-8/'+str(t)+'tdyn/'

    if not os.path.exists(SAVE_PATH): os.makedirs(SAVE_PATH)
    s.save(SAVE_PATH + plot[0].replace(".","_") + '_' + str(plot[1]+t).zfill(4) + '.pdf')
    
    