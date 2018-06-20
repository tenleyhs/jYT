import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
execfile('/pfs/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

"""
fig, ax = plt.subplots()
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p10_b3_0_0040_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C0', ls='--', label='128, chk40')
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p10_b3_0_0060_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C0', label='128, chk60')
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p10_b3_0_256_0040_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C1', ls='--', label='256, chk40')
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p10_b3_0_256_0057_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C1', label='256, chk57')
ax.set_xlim(-2, 2)
ax.legend()
fig.tight_layout()
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/temp/compare_resolution_p10_b3_0.pdf')
plt.close('all')

fig, ax = plt.subplots()
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p10_b2_0_0040_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C0', ls='--', label='128, chk40')
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p10_b2_0_0075_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C0', label='128, chk75')
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p10_b2_0_256_0040_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C1', ls='--', label='256, chk40')
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p10_b2_0_256_0080_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C1', label='256, chk80')
ax.set_xlim(-2, 2)
ax.legend()
fig.tight_layout()
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/temp/compare_resolution_p10_b2_0.pdf')
"""

fig, ax = plt.subplots()
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p16_b3_0_0040_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C0', ls='--', label='24k, chk40')
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p16_b3_0_0060_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C0', label='24k, chk60')
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p16_b3_0_48k_0035_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C1', ls='--', label='48k, chk35')
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p16_b3_0_48k_0040_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C1', ls=':', label='48k, chk40')
x, y = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/results/dmdts/data/m1_0_p16_b3_0_48k_0048_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C1', ls=':', label='48k, chk48')
ax.set_xlim(-2, 2)
ax.legend()
fig.tight_layout()
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/temp/compare_resolution_p16_b3_0.pdf')

"""
import yt
import colormaps as cmaps

data = [
    ['/pfs/lawsmith/FLASH4.3/runs/m1.0_p10_b2.0/multitidal_hdf5_plt_cnt_0400', 'm1_0_p10_b2_0_0400'],
    ['/pfs/lawsmith/FLASH4.3/runs/m1.0_p10_b2.0_256/multitidal_hdf5_plt_cnt_0400', 'm1_0_p10_b2_0_256_0400'],
    ['/pfs/lawsmith/FLASH4.3/runs/m1.0_p10_b3.0/multitidal_hdf5_plt_cnt_0400', 'm1_0_p10_3_0_0400'],
    ['/pfs/lawsmith/FLASH4.3/runs/m1.0_p10_b3.0_256/multitidal_hdf5_plt_cnt_0400', 'm1_0_p10_b3_0_256_0400']
]

zlim = (1e-4, 1e0)
width = (50, 'rsun')
var = 'dens'
log_tf = True

for dat in data:
    ds = yt.load(dat[0])
    s = yt.SlicePlot(ds, 'z', var, width=width)
    s.set_log(var, log_tf)
    s.annotate_contour("density", ncont=3, clim=(1e-3,1e0), label=True,
                            plot_args={"colors":"white","linewidths":1})
    s.set_zlim(var, zlim[0], zlim[1])
    s.set_cmap(var, cmaps.viridis)
    s.annotate_scale(unit='rsun', text_args={"size":48})
    s.hide_axes()
    s.hide_colorbar()
    s.set_figure_size(12)
    s.save('/pfs/lawsmith/FLASH4.3/runs/results/temp/' + dat[1] + '_' + str(width[0]) + 'rsun' + '.pdf')
"""
