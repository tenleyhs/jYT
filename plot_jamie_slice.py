"""
simple slices
"""
import yt
import colormaps as cmaps

plots = [
    ['m1.0_p1_b2.0', '0023', 0.8*80.78, 1e-5*80.78, 'age = 0Gyr'],
    ['m1.0_p10_b2.0', '0023', 0.8*169.88, 1e-5*169.88, 'age = 4.8Gyr'],
    ['m1.0_p16_b2.0', '0023', 0.8*500, 1e-5*500, 'age = 8.4Gyr'],
]

var = 'dens'
width = (10, 'rsun')
log_tf = True
set_zlim_tf = True
contour_tf = True
sz = 48

LOAD_FILES = []
for plot in plots:
    ds = yt.load('/storage/dark/lawsmith/FLASH4.3/runs/'+plot[0]+'/multitidal_hdf5_chk_'+plot[1])
    s = yt.SlicePlot(ds, 'z', var, width=width)
    s.set_log(var, log_tf)
    if contour_tf:
        s.annotate_contour("density", ncont=3, clim=(1e-2,1e0), label=False,
                            plot_args={"colors":"white","linewidths":2,"linestyles":"-"},
                            text_args={"size":48})
    if set_zlim_tf:
        s.set_zlim(var, plot[3], plot[2])
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
    s.save('/groups/dark/lawsmith/results/paper/' + plot[0].replace(".","_") + '_' + plot[1] + '.pdf')
