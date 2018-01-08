'''
makes a bunch of slices in parallel in prep for a movie
then can run:
ffmpeg -r 24 -f image2 -i multitidal_hdf5_plt_cnt_%04d_Slice_z_dens.png -b:v 4M movie.mpg
or can use quicktime 7 like James suggested
'''

dir = 'm1.0_p10_b2.0_256'
var = 'dens'
log_tf = True
width = (10, 'rsun')
zmax = 169.88  # p1:80.78, p10:169.88, p16:500   # max density of this profile
zmin = 1e-5 * zmax
set_zlim_tf = True

#LOAD_FILES = '/nobackup/jlawsmit/' + dir + '/multitidal_hdf5_plt_cnt_*'
#LOAD_FILES = '/nobackup/jlawsmit/' + dir + '/multitidal_hdf5_plt_cnt_[0-1][2-9][0-9][0-9]'
#SAVE_PATH = '/nobackup/jlawsmit/' + dir + '/plots/'

#LOAD_FILES = '/pfs/lawsmith/FLASH4.3/runs/' + dir + '/multitidal_hdf5_plt_cnt_*'
LOAD_FILES = '/pfs/lawsmith/FLASH4.3/runs/' + dir + '/multitidal_hdf5_plt_cnt_0[2-9][0-9][0-9]'
SAVE_PATH = '/pfs/lawsmith/FLASH4.3/runs/' + dir + '/plots/'

import yt
import colormaps as cmaps
yt.enable_parallelism()

ts = yt.DatasetSeries(LOAD_FILES)

for ds in ts.piter():
	s = yt.SlicePlot(ds, 'z', var, width=width)
	s.set_log(var, log_tf)
	if set_zlim_tf:
		s.set_zlim(var, zmin, zmax)
	s.set_cmap(var, cmaps.viridis)
	s.annotate_timestamp(time_unit='s')
	s.annotate_scale(unit='rsun')
	#s.set_axes_unit('rsun')
	#s.set_font({'size':24})
	s.hide_axes()
	s.hide_colorbar()
	s.set_figure_size(12)
	#s.set_buff_size(2000)
	s.save(SAVE_PATH)
