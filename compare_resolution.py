import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
#execfile('/pfs/lawsmith/jYT/my_settings.py')
execfile('/groups/dark/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

## TODO can't I do this locally?

"""
ds = [
    ['m1.0_p10_b3.0',		'0043',     20,     -2,	2,        0,  0,  0],
    ['m1.0_p10_b3.0',		'0030',     1,     -2,	2,        0,  0,  0],
    ['43_b2.0_6k',		'0040',     1,     -2,	2,        0,  0,  0],
    ['43_b2.0_6k',		'0060',     1,     -2,	2,        0,  0,  0],
    ['43_b2.0_6k',		'0080',     1,     -2,	2,        0,  0,  0],
    ['43_b2.0_12k',		'0041',     1,     -2,	2,        0,  0,  0]
    ]

fig, ax = plt.subplots()
x, y = np.loadtxt('/pfs/lawsmith/results/dmdts/data/m1_0_p10_b2_0_0075_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C0', label='24k, chk75')
x, y = np.loadtxt('/pfs/lawsmith/results/dmdts/data/duplicate_m1_0_p10_b2_0_0075_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C1', ls='--', label='24k, chk75 (duplicate)')
ax.set_xlim(-2, 2)
ax.legend()
fig.tight_layout(pad=0.3)
fig.savefig('/pfs/lawsmith/results/temp/compare_resolution_p10_b2_0_duplicate.pdf')
"""
"""
fig, ax = plt.subplots()
x, y = np.loadtxt('/groups/dark/lawsmith/results/m1_0_p10_b3_0_0043_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C1', ls='--', label='m1.0_p10_b3.0, chk43')
ax.set_xlim(-2, 2)
ax.legend()
fig.tight_layout(pad=0.3)
fig.savefig('/groups/dark/lawsmith/results/compare_resolution/compare_resolution_m1_0_p10_b3_0.pdf')
"""

fig, ax = plt.subplots()
x, y = np.loadtxt('/groups/dark/lawsmith/results/43_b2_0_12k_0041_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C1', ls='--', label='4/3, 12k, chk41')
x, y = np.loadtxt('/groups/dark/lawsmith/results/43_b2_0_12k_0080_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C2', ls='--', label='4/3, 12k, chk80')
x, y = np.loadtxt('/groups/dark/lawsmith/Guillochon2013_dmdts/4-3/2.000.dat', unpack=True)
ax.plot(x, y, c='C3', label='4/3, GRR2013')
ax.set_xlim(-2, 2)
ax.legend()
fig.tight_layout(pad=0.3)
fig.savefig('/groups/dark/lawsmith/results/compare_resolution/compare_resolution_43_b2_0_12k.pdf')
"""
fig, ax = plt.subplots()
x, y = np.loadtxt('/pfs/lawsmith/results/dmdts/data/43_b2_0_24k_0080_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C0', ls='--', label='24k, chk80')
x, y = np.loadtxt('/pfs/lawsmith/results/dmdts/data/43_b2_0_48k_0080_ev.dat', skiprows=1, unpack=True)
ax.plot(x, y, c='C1', ls='--', label='48k, chk80')
x, y = np.loadtxt('/pfs/lawsmith/dmdts/4-3/2.000.dat', unpack=True)
ax.plot(x, y, c='C2', label='4/3, GRR2013')
ax.set_xlim(-2, 2)
ax.legend()
fig.tight_layout(pad=0.3)
fig.savefig('/pfs/lawsmith/results/temp/compare_resolution_43_b2_0.pdf')
"""
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
fig.tight_layout(pad=0.3)
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
fig.tight_layout(pad=0.3)
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/temp/compare_resolution_p10_b2_0.pdf')
"""
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
fig.tight_layout(pad=0.3)
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/temp/compare_resolution_p16_b3_0.pdf')
"""
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
