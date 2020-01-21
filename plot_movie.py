'''
run like:
mpirun -np 8 python ~/jYT/plot_movie.py --cluster=fend --run=m1.0_p16_b3.0 --var=dens --width=1000
mpirun -np 8 python ~/jYT/plot_movie.py --run=1e3_m1.0_p16_b5.0_48k --width=200 --chkplt=chk --denscut=True --var='he4 ' --log=False
or preferably with:
the runfile in this directory, in e.g. the plot1/ run.

makes a bunch of slices in parallel in prep for a movie
then can run:
ffmpeg -f image2 -framerate 24 -i multitidal_hdf5_chk_%04d_Slice_z_dens.png -b:v 4M movie.mpg
ffmpeg -f image2 -framerate 24 -i multitidal_hdf5_plt_cnt_%04d_Slice_z_dens.png -b:v 4M movie.mpg
ffmpeg -pattern_type glob -i "multitidal_hdf5_plt_cnt_*_Slice_z_dens.png" -b:v 4M movie.mpg
'''
import yt
import os
import argparse, sys

parser = argparse.ArgumentParser()
parser.add_argument('--cluster', help='cluster, e.g., hyades/fend/pfe', type=str, default='fend')
parser.add_argument('--run', help='run to plot, e.g., m1.0_p16_b3.0', type=str)
parser.add_argument('--var', help='variable to plot, e.g., dens or he4 ', type=str, default='dens')
parser.add_argument('--width', help='width of plot in rsun, e.g., 1000', type=int, default=99999)
parser.add_argument('--zlim', help='whether to impose zlim, e.g. True/False', type=bool, default=False)
parser.add_argument('--zmax', help='maximum of variable to plot, e.g., 80.78/169.88/500', type=float, default=500)
parser.add_argument('--zmin', help='minimum of variable to plot, e.g., 80.78/169.88/500', type=float, default=500)
parser.add_argument('--factor', help='sets minimum of variable to plot', type=float, default=1e-5)
parser.add_argument('--chkplt', help='plot from chks or plts', type=str, default='chk')
parser.add_argument('--vel', help='annotate velocity?', type=bool, default=False)
parser.add_argument('--denscut', help='do density cut?', type=bool, default=False)
parser.add_argument('--denscutval', help='value of density cut', type=float, default=1e-5)
parser.add_argument('--log', help='log scale?', type=bool, default=True)
parser.add_argument('--files', help='which files to plot, e.g. [0-9][0-9][0-9]0 or *', type=str, default='*')
parser.add_argument('--proj', help='projection plot?', type=bool, default=False)
args = parser.parse_args()

if args.cluster == 'hyades':
	clusterdir = '/pfs/lawsmith/FLASH4.3/runs/'
	savepath = '/pfs/lawsmith/hyades_slices/'
elif args.cluster == 'fend':
	clusterdir = '/storage/dark/lawsmith/FLASH4.3/runs/'
	savepath = '/groups/dark/lawsmith/fend_slices/'
	#clusterdir = '/storage/dark/lawsmith/FLASH4.5/runs_bondi/'
	#savepath = '/groups/dark/lawsmith/bondi_slices/'
elif args.cluster == 'pfe':
	clusterdir = '/nobackup/jlawsmit/'
	savepath = '/nobackup/jlawsmit/pfe_slices/'
elif args.cluster == 'lux':
	clusterdir = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/'
	savepath = '/data/groups/ramirez-ruiz/lawsmith/FLASH4.3/lux_slices/'
else:
	print("haven't set up this cluster in script yet")
	exit()

#if not os.path.exists(savepath + args.run + '/' + args.var + '_' + str(args.width) + 'rsun/'):
#	os.makedirs(savepath + args.run + '/' + args.var + '_' + str(args.width) + 'rsun/')

if args.chkplt == 'chk':
	LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_chk_' + args.files
elif args.chkplt == 'plt':
	LOAD_FILES = clusterdir + args.run + '/multitidal_hdf5_plt_cnt_' + args.files
#LOAD_FILES = clusterdir + args.run + '/bondi_hdf5_chk_*'
#LOAD_FILES = clusterdir + args.run + '/bondi_hdf5_plt_cnt_*'
yt.enable_parallelism()
ts = yt.DatasetSeries(LOAD_FILES)

for ds in ts.piter():
	if args.denscut:
		ad = ds.all_data()
		dense_ad = ad.cut_region(['obj["dens"] > ' + str(args.denscutval) ]) #1e-5'])
		if args.width == 99999:
			s = yt.SlicePlot(ds, 'z', args.var, data_source=dense_ad)
		else:
			s = yt.SlicePlot(ds, 'z', args.var, data_source=dense_ad, width=(args.width, 'rsun'))
		# for some reason the log swtich doesn't seem to work here. have to hard code it
		s.set_log(args.var, False)
	else:
		if args.width == 99999:
			if args.proj:
				s = yt.ProjectionPlot(ds, 'z', args.var)
			else:
				s = yt.SlicePlot(ds, 'z', args.var)

		else:
			s = yt.SlicePlot(ds, 'z', args.var, width=(args.width, 'rsun'))
		s.set_log(args.var, args.log)

	if args.vel:
		s.annotate_velocity()

	if args.zlim:
		#s.set_zlim(args.var, args.factor * args.zmax, args.zmax)
		s.set_zlim(args.var, args.zmin, args.zmax)


	#s.set_log(args.var, False)    # TODO jamie temp
	#s.set_cmap(var, cmaps.viridis)
	s.annotate_timestamp(time_unit='day')
	s.annotate_scale(unit='rsun')
	#s.set_axes_unit('rsun')
	#s.set_font({'size':24})
	#s.hide_axes()
	#s.hide_colorbar()
	#s.set_figure_size(12)
	#s.set_buff_size(2000)

	#s.annotate_line_integral_convolution('velx', 'vely')#, lim=(0.5,0.65))
	if args.zlim:
		if args.denscut:
			s.save(savepath +args.run+'/'+args.var+'_'+str(args.width)+'rsun'+'_'+str(args.zmax)+'_'+str(args.factor)+'_denscut'+str(args.denscutval)+'/')
		else:
			s.save(savepath +args.run+'/'+args.var+'_'+str(args.width)+'rsun'+'_'+str(args.zmax)+'_'+str(args.factor)+'/')
	else:
		if args.denscut:
			s.save(savepath +args.run+'/'+args.var+'_'+str(args.width)+'rsun'+'_nozlim'+'_denscut'+str(args.denscutval)+'/')
		else:
			if args.proj:
				s.save(savepath +args.run+'/'+'proj_'+args.var+'_'+str(args.width)+'rsun'+'_nozlim/')
			else:
				s.save(savepath +args.run+'/'+args.var+'_'+str(args.width)+'rsun'+'_nozlim/')
