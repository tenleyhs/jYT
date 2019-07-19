"""
plot elements in slices, following Enrico meeting.
"""

#dir = 'm3.0_p16_b4.0_48k'
#width = (30, 'rsun')
dir = 'm1.0_p16_b3.0_96k_1000'
width = (10, 'rsun')
zmax = 1
zmin = 0

var = 'he4 '
dens_cut_tf = True

#var = 'c12 '
#dens_cut_tf = True

#var = 'n14 '
#dens_cut_tf = False

LOAD_FILES = '/storage/dark/lawsmith/FLASH4.3/runs/' + dir + '/multitidal_hdf5_chk_0021'
#LOAD_FILES = '/storage/dark/lawsmith/FLASH4.3/runs/' + dir + '/multitidal_hdf5_chk_0000'
SAVE_PATH = '/groups/dark/lawsmith/results/paper/' + dir + '/' + var[:3]+ '/'

import yt
import colormaps as cmaps
yt.enable_parallelism()

ts = yt.DatasetSeries(LOAD_FILES)

for ds in ts.piter():
	ad = ds.all_data()
	if dens_cut_tf:
		dense_ad = ad.cut_region(['obj["dens"] > 1e-4'])
		s = yt.SlicePlot(ds, 'z', var, data_source=dense_ad, width=width)
	else:
		s = yt.SlicePlot(ds, 'z', var, width=width)
	s.set_log(var, False)
	#s.set_zlim(var, zmin, zmax)
	s.set_cmap(var, cmaps.viridis)
	#s.annotate_timestamp(time_unit='s')
	##s.annotate_scale(unit='rsun', text_args={"size":48})
	#s.set_axes_unit('rsun')
	#s.set_font({'size':24})
	s.hide_axes()
	s.hide_colorbar()
	s.set_figure_size(12)
	#s.set_buff_size(8000)
	s.save(SAVE_PATH)
