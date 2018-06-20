'''
mpirun -np 8 python movie.py
makes a bunch of slices in parallel in prep for a movie
then can run:
ffmpeg -r 24 -f image2 -i multitidal_hdf5_plt_cnt_%04d_Slice_z_dens.png -b:v 4M movie.mpg
or can use quicktime 7 like James suggested
'''

cluster = '/groups/dark/lawsmith/FLASH4.3_copy/runs/'
#cluster = '/pfs/lawsmith/FLASH4.3/runs/'

dir = 'm1.0_p16_b3.0'

var = 'pres'
log_tf = True
width = (50, 'rsun')
set_zlim_tf = False
zmax = 169.88   # p1:80.78, p10:169.88, p16:500
zmin = 1e-5 * zmax


LOAD_FILES = cluster + dir + '/multitidal_hdf5_chk_*'
SAVE_PATH = cluster + dir + '/' + var + '_' + str(width[0]) + 'rsun/'

import yt
#import colormaps as cmaps
yt.enable_parallelism()

ts = yt.DatasetSeries(LOAD_FILES)

for ds in ts.piter():
	s = yt.SlicePlot(ds, 'z', var, width=width)
	s.set_log(var, log_tf)
	if set_zlim_tf: s.set_zlim(var, zmin, zmax)
	#s.set_cmap(var, cmaps.viridis)
	s.annotate_timestamp(time_unit='s')
	s.annotate_scale(unit='rsun')
	#s.set_axes_unit('rsun')
	#s.set_font({'size':24})
	#s.hide_axes()
	#s.hide_colorbar()
	#s.set_figure_size(12)
	#s.set_buff_size(2000)
	#s.annotate_velocity()
	#s.annotate_line_integral_convolution('velx', 'vely')#, lim=(0.5,0.65))
	s.save(SAVE_PATH)
