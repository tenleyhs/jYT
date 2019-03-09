'''
run like:
mpirun -np 8 python jYT/plot_movie.py --cluster=fend --run=m1.0_p16_b3.0 --var=dens --width=1000

makes a bunch of slices in parallel in prep for a movie
then can run:
ffmpeg -r 24 -i multitidal_hdf5_chk_%04d_Slice_z_dens.png -b:v 4M movie.mp4
ffmpeg -r 24 -i multitidal_hdf5_plt_cnt_0%02d0_Slice_z_dens.png -b:v 4M movie.mp4
or can use quicktime 7 like James suggested
'''
import yt
import os
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='hyades')
parser.add_argument('--run', help='run to plot, e.g., m1.0_p16_b3.0', type=str)
parser.add_argument('--var', help='variable to plot, e.g., dens', type=str, default='dens')
parser.add_argument('--width', help='width of plot in rsun, e.g., 1000', type=int, default=1000)
parser.add_argument('--zmax', help='maximum of variable to plot, e.g., 80.78/169.88/500', type=float, default=500)
args = parser.parse_args()

if args.cluster == 'hyades':
	clusterdir = '/pfs/lawsmith/FLASH4.3/runs/'
	savepath = '/pfs/lawsmith/hyades_slices/'
elif args.cluster == 'fend':
	clusterdir = '/storage/dark/lawsmith/FLASH4.3/runs/'
	savepath = '/groups/dark/lawsmith/fend_slices/'
elif args.cluster == 'pfe':
	clusterdir = '/nobackup/jlawsmit/'
	savepath = '/nobackup/jlawsmit/pfe_slices/'
else:
	print("haven't set up this cluster in script yet")
	exit()

#if not os.path.exists(savepath + args.run + '/' + args.var + '_' + str(args.width) + 'rsun/'):
#	os.makedirs(savepath + args.run + '/' + args.var + '_' + str(args.width) + 'rsun/')

LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_chk_*'
yt.enable_parallelism()
ts = yt.DatasetSeries(LOAD_FILES)

for ds in ts.piter():
	s = yt.SlicePlot(ds, 'z', args.var, width=(args.width, 'rsun'))
	s.set_log(args.var, True)
	#s.set_zlim(args.var, 1e-5 * args.zmax, args.zmax)
	#s.set_cmap(var, cmaps.viridis)
	s.annotate_timestamp(time_unit='s')
	s.annotate_scale(unit='rsun')
	#s.set_axes_unit('rsun')
	#s.set_font({'size':24})
	#s.hide_axes()
	#s.hide_colorbar()
	#s.set_figure_size(12)
	#s.set_buff_size(2000)
	s.annotate_velocity()
	#s.annotate_line_integral_convolution('velx', 'vely')#, lim=(0.5,0.65))
	s.save(savepath + args.run + '/' + args.var + '_' + str(args.width) + 'rsun/')
