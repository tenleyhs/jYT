"""
plot elements in slices, following Enrico meeting.
"""

dir = 'm1.0_p16_b2.0'
var = 'he4 '
width = (10, 'rsun')
log_tf = False
set_zlim_tf = True
zmax = 1
zmin = 0

#LOAD_FILES = '/pfs/lawsmith/FLASH4.3/runs/' + dir + '/multitidal_hdf5_plt_cnt_0[1-2][0-9]0'
LOAD_FILES = '/pfs/lawsmith/FLASH4.3/runs/' + dir + '/multitidal_hdf5_plt_cnt_0240'
SAVE_PATH = '/pfs/lawsmith/FLASH4.3/runs/' + dir + '/he4/'

import yt
import colormaps as cmaps
yt.enable_parallelism()

ts = yt.DatasetSeries(LOAD_FILES)

for ds in ts.piter():
	ad = ds.all_data()
	dense_ad = ad.cut_region(['obj["dens"] > 1e-4'])
	s = yt.SlicePlot(ds, 'z', var, data_source=dense_ad, width=width)
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
