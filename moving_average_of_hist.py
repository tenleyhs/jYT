import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
execfile('/pfs/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

M_bh = 1e6*M_sun

e, dm = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/m1.0_p1_b1.0_pfe/b10000_ev_bhbound_histogram_multitidal_hdf5_chk_0100.dat', skiprows=4)
de = e[1]-e[0]
dm_de = dm/de
log_dm_de = np.log10(dm_de)
e_bound = e[np.where(e<0.)]
t = 2.*np.pi*G*M_bh/((2*np.abs(e_bound))**(3./2))
de_dt = (1./3)*((2*np.pi*G*M_bh)**(2./3))*t**(-5./3)
dm_de_bound = dm_de[np.where(e<0)]
mdot = dm_de_bound*de_dt

log_t_yr = np.log10(t/yr)
mdot_moyr = mdot*yr/M_sun
log_mdot_moyr = np.log10(mdot_moyr)

# plot unsmoothed
fig, ax = plt.subplots()
ax.scatter(log_t_yr, log_mdot_moyr, c='C0', rasterized=True,
    alpha=0.5, s=1, edgecolors='none')

x_arr = []
y_arr = []
for x in np.linspace(-2, 1.95, num=80):
    ixs = np.logical_and(log_t_yr >= x, log_t_yr < x+0.05)
    y = np.mean(mdot_moyr[ixs])
    x_arr.append(x+0.025)
    y_arr.append(np.log10(y))

ax.plot(x_arr, y_arr, c='C1')
y_arr = np.array(y_arr)
print max(y_arr[np.isfinite(y_arr)])

ax.set_xlim(-2, 2)
ax.set_ylabel(r'$\log\ \dot M\ {\rm [M_\odot/yr]}$')
ax.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
fig.tight_layout(pad=0.3)
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/temp/mdot_p1_b1_0_pfe.pdf')
