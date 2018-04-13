"""
simple slices
"""
import yt
import colormaps as cmaps

plots = [
    ['m1.0_p1_b2.0', '0230', 0.8*80.78, 1e-5*80.78],
    ['m1.0_p10_b2.0_256', '0230', 0.8*169.88, 1e-5*169.88],
    ['m1.0_p16_b2.0', '0230', 0.8*500, 1e-5*500],
]

var = 'dens'
width = (10, 'rsun')
log_tf = True
set_zlim_tf = True

LOAD_FILES = []
for plot in plots:
    ds = yt.load('/pfs/lawsmith/FLASH4.3/runs/'+plot[0]+'/multitidal_hdf5_plt_cnt_'+plot[1])
    s = yt.SlicePlot(ds, 'z', var, width=width)
    s.set_log(var, log_tf)
    if set_zlim_tf:
        s.set_zlim(var, plot[3], plot[2])
    s.set_cmap(var, cmaps.viridis)
    s.annotate_timestamp(time_unit='s')
    s.annotate_scale(unit='rsun')
    #s.set_axes_unit('rsun')
    #s.set_font({'size':24})
    s.hide_axes()
    s.hide_colorbar()
    s.set_figure_size(12)
    #s.set_buff_size(2000)
    s.save('/pfs/lawsmith/FLASH4.3/runs/results/' + plot[0].replace(".","_") + '_' + plot[1] + '.pdf')
