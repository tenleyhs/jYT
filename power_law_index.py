import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from astropy.io import ascii
from scipy.interpolate import UnivariateSpline
#execfile('/nobackup/jlawsmit/jYT/my_settings.py')
execfile('/pfs/lawsmith/jYT/my_settings.py')
plt.rcParams['legend.fontsize'] = 16
plt.rcParams['font.size'] = 18

# PARAMS
M_bh = 1e6*M_sun
lw = 1.5
s = 20

# CHANGE EACH TIME
fname = 'power_law_index_inf'

ages = [
        [
    ['m1.0_p1_b1.0',1.0,'multitidal_hdf5_chk_0100',300,0.],
    ['m1.0_p1_b2.0',2.0,'multitidal_hdf5_chk_0080',250,-0.1],
    ['m1.0_p1_b3.0',3.0,'multitidal_hdf5_chk_0060',250,0.3]
        ],
        [
    ['m1.0_p10_b1.0',1.0,'multitidal_hdf5_chk_0100',400,-0.4],
    ['m1.0_p10_b1.5',1.5,'multitidal_hdf5_chk_0090',200,-0.15],
    ['m1.0_p10_b2.0_128',2.0,'multitidal_hdf5_chk_0075',200,0.25],
    #['m1.0_p10_b2.5',2.5,'multitidal_hdf5_chk_0080'],
    ['m1.0_p10_b3.0',3.0,'multitidal_hdf5_chk_0060',200,0.4],
    ['m1.0_p10_b4.0',4.0,'multitidal_hdf5_chk_0050',200,0.25],
    #['m1.0_p10_b5.0',5.0,'multitidal_hdf5_chk_0060']
        ],
        [
    ['m1.0_p16_b2.0',2.0,'multitidal_hdf5_chk_0075',300,3],
    ['m1.0_p16_b3.0',3.0,'multitidal_hdf5_chk_0060',500,0],
    ['m1.0_p16_b4.0',4.0,'multitidal_hdf5_chk_0050',300,0.1],
    ['m1.0_p16_b5.0',5.0,'multitidal_hdf5_chk_0050',300,0.1]
        ]
       ]


labels = ['age=0Gyr', 'age=4.8Gyr', 'age=8.4Gyr']

# limits by eye for t in dmdt plot. from jupyter
min_log_t = 6.5
max_log_t = 9.0

# for dmde's
def smooth(x,window_len=11,window='hanning'):
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."
	if window_len<3:
		return x
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
	s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
	if window == 'flat': #moving average
		w=np.ones(window_len,'d')
	else:
		w=eval('np.'+window+'(window_len)')
	y=np.convolve(w/w.sum(),s,mode='same')
	return y[window_len:-window_len+1]


fig, ax = plt.subplots()

for i, age in enumerate(ages):
    b_array = []
    n_array = []
    for p in age:
        e, dm = np.loadtxt('/pfs/lawsmith/FLASH4.3/runs/' + p[0] + '/b10000_ev_bhbound_histogram_' + p[2] + '.dat', skiprows=4)
        de = e[1]-e[0]
        dm_de = dm/de

        # calculate mdot
        e_bound = e[np.where(e<0.)]
        t = 2.*np.pi*G*M_bh/((2*np.abs(e_bound))**(3./2))
        de_dt = (1./3)*((2*np.pi*G*M_bh)**(2./3))*t**(-5./3)
        dm_de_bound = dm_de[np.where(e<0)]
        mdot = dm_de_bound*de_dt

        # take logs. better to smooth what I'm actually plotting (log-linear or log-log)
        log_t = np.log10(t)
        log_t_yr = np.log10(t/yr)
        log_mdot = np.log10(mdot)
        log_mdot_moyr = np.log10(mdot*yr/M_sun)
        log_dm_de = np.log10(dm_de)

        slog_mdot_moyr = smooth(log_mdot_moyr, window_len=p[3])
        sel = np.where(log_t_yr<p[4])[0]
        x = log_t_yr[sel]
        y = slog_mdot_moyr[sel]

        #TODO better slope function
        ninf = (y[-2] - y[-1]) / (x[-2] - x[-1])
        ninf = (y[-2] - y[-1]) / (x[-2] - x[-1])

        b_array.append(p[1])
        n_array.append(ninf)

    ax.plot(b_array, n_array, lw=lw, label=labels[i], alpha=0.5)
    ax.scatter(b_array, n_array, s=s)


#ax.set_xlim(1.0, 4.0)
#ax.set_ylim(-2, 1)
ax.set_ylabel(r'$n_\infty$')
ax.set_xlabel(r'$\log\ t\ \mathrm{[yr]}$')
ax.legend()
fig.tight_layout()
fig.savefig('/pfs/lawsmith/FLASH4.3/runs/results/'+fname+'.png')
